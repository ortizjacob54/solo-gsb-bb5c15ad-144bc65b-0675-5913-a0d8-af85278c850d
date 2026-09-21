"""FastAPI 应用入口。

启动时（lifespan）自动完成数据库初始化；接口：
    POST /digests        提交序列与酶列表，执行酶切并保存，返回完整结果
    GET  /digests        列出全部提交的摘要
    GET  /digests/{id}   按 id 查询保存的完整酶切结果
    GET  /health         健康检查
"""

from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException

from .db import Database
from .digest import run_digest
from .schemas import DigestRequest, DigestResponse, DigestSummary

db: Database


@asynccontextmanager
async def lifespan(app: FastAPI):
    global db
    db = Database(os.environ.get("DATABASE_PATH", "data/digest.db"))
    yield


app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description="DNA 限制性酶切模拟：支持线性/环形序列、IUPAC 简并碱基、"
                "双链切口偏移、黏端/平端分析与片段生成。",
    lifespan=lifespan,
)


def _build_response(record: dict) -> dict:
    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence": record["sequence"],
        "sequence_length": len(record["sequence"]),
        "enzymes": record["enzymes"],
        "cuts": record["result"]["cuts"],
        "fragments": record["result"]["fragments"],
        "created_at": record["created_at"],
    }


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/digests", status_code=201, response_model=DigestResponse)
def create_digest(request: DigestRequest) -> dict:
    # Pydantic 已在此前完成全部校验（含切口偏移合法性），
    # 非法请求在写库之前即以 422 拒绝，保证整次保存被拒绝。
    result = run_digest(request.sequence, request.topology, request.enzymes)
    submission_id, created_at = db.save_submission(
        request.sequence, request.topology, request.enzymes, result
    )
    return _build_response(
        {
            "id": submission_id,
            "sequence": request.sequence,
            "topology": request.topology,
            "created_at": created_at,
            "enzymes": [e.model_dump() for e in request.enzymes],
            "result": result,
        }
    )


@app.get("/digests", response_model=list[DigestSummary])
def list_digests() -> list[dict]:
    return db.list_submissions()


@app.get("/digests/{submission_id}", response_model=DigestResponse)
def get_digest(submission_id: int) -> dict:
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")
    return _build_response(record)
