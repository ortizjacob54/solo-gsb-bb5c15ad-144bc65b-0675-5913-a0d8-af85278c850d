"""变异影响分析 API 测试：新增/消失/末端改变的切口、片段长度差异与校验。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
# 简并模式 GANNTC：识别位点中部两个碱基任意，
# 便于在不破坏位点的前提下改变黏端突出序列
ND_ENZ = {"name": "Nd", "pattern": "GANNTC", "cut_top": 1, "cut_bottom": 5}


def _create(client, **overrides) -> int:
    payload = {
        "sequence": "AAAAGAATTCTTTT",  # EcoRI 位点起于 4，切口 (5, 9)
        "topology": "linear",
        "enzymes": [ECORI],
    }
    payload.update(overrides)
    resp = client.post("/digests", json=payload)
    assert resp.status_code == 201
    return resp.json()["id"]


def _analyze(client, sub_id, variants):
    return client.post(
        f"/digests/{sub_id}/variant-analysis", json={"variants": variants}
    )


def test_variant_creates_new_cut(client):
    # 原序列 GAATTT 不含 EcoRI 位点；position 9 T->C 后变为 GAATTC
    sub_id = _create(client, sequence="AAAAGAATTTTTTT")
    resp = _analyze(client, sub_id, [{"position": 9, "base": "C"}])
    assert resp.status_code == 200
    body = resp.json()
    assert body["submission_id"] == sub_id
    assert body["topology"] == "linear"
    assert body["sequence_length"] == 14
    assert body["mutated_sequence"] == "AAAAGAATTCTTTT"
    assert body["variants_applied"] == [{"position": 9, "ref_base": "T", "base": "C"}]

    assert len(body["added_cuts"]) == 1
    cut = body["added_cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (5, 9)
    assert cut["end_type"] == "5prime"
    assert cut["overhang"] == "AATT"
    assert cut["enzymes"] == ["EcoRI"]
    assert body["removed_cuts"] == []
    assert body["changed_cuts"] == []

    assert body["fragment_lengths"] == {
        "baseline_lengths": [14],
        "variant_lengths": [5, 9],
        "added_lengths": [5, 9],
        "removed_lengths": [14],
    }


def test_variant_removes_existing_cut(client):
    sub_id = _create(client)  # GAATTC 起于 4
    resp = _analyze(client, sub_id, [{"position": 9, "base": "T"}])
    assert resp.status_code == 200
    body = resp.json()
    assert body["added_cuts"] == []
    assert len(body["removed_cuts"]) == 1
    assert (body["removed_cuts"][0]["top_cut"], body["removed_cuts"][0]["bottom_cut"]) == (5, 9)
    assert body["changed_cuts"] == []
    assert body["fragment_lengths"]["baseline_lengths"] == [5, 9]
    assert body["fragment_lengths"]["variant_lengths"] == [14]
    assert body["fragment_lengths"]["added_lengths"] == [14]
    assert body["fragment_lengths"]["removed_lengths"] == [5, 9]


def test_variant_changes_overhang_only(client):
    # GANNTC 中 N 匹配任意碱基：position 6 A->C 不破坏位点，但突出端 AATT -> ACTT
    sub_id = _create(client, enzymes=[ND_ENZ])
    resp = _analyze(client, sub_id, [{"position": 6, "base": "C"}])
    assert resp.status_code == 200
    body = resp.json()
    assert body["added_cuts"] == [] and body["removed_cuts"] == []
    assert len(body["changed_cuts"]) == 1
    change = body["changed_cuts"][0]
    assert (change["top_cut"], change["bottom_cut"]) == (5, 9)
    assert change["changes"] == ["overhang"]
    assert change["before"] == {"end_type": "5prime", "overhang": "AATT"}
    assert change["after"] == {"end_type": "5prime", "overhang": "ACTT"}
    # 切口坐标未变，片段长度也不变
    assert body["fragment_lengths"]["added_lengths"] == []
    assert body["fragment_lengths"]["removed_lengths"] == []


def test_multiple_variants_applied_together(client):
    # 位点 1（GAATTC@4）被破坏，同时在 14 处新建一个位点
    sub_id = _create(client, sequence="AAAAGAATTCTTTTGAATTT")
    # 请求中故意乱序传入，响应里的 variants_applied 应按位置排序
    resp = _analyze(client, sub_id, [{"position": 19, "base": "C"},
                                     {"position": 9, "base": "T"}])
    assert resp.status_code == 200
    body = resp.json()
    assert [v["position"] for v in body["variants_applied"]] == [9, 19]
    assert body["variants_applied"][0]["ref_base"] == "C"
    assert body["variants_applied"][1]["ref_base"] == "T"
    assert [c["top_cut"] for c in body["removed_cuts"]] == [5]
    assert [c["top_cut"] for c in body["added_cuts"]] == [15]
    # 长度多重集相同（5/15 互换），差异为空
    assert body["fragment_lengths"]["baseline_lengths"] == [5, 15]
    assert body["fragment_lengths"]["variant_lengths"] == [15, 5]
    assert body["fragment_lengths"]["added_lengths"] == []
    assert body["fragment_lengths"]["removed_lengths"] == []


def test_circular_variant_removes_origin_spanning_cut(client):
    # GAATTC 跨原点（start=8，占位置 8,9,0,1,2,3）；position 0 A->G 破坏位点
    sub_id = _create(client, sequence="ATTCGGGGGA", topology="circular")
    resp = _analyze(client, sub_id, [{"position": 0, "base": "G"}])
    assert resp.status_code == 200
    body = resp.json()
    assert [c["top_cut"] for c in body["removed_cuts"]] == [9]
    assert body["added_cuts"] == []
    # 单切口开环与无切口均为一条全长片段，长度差为空
    assert body["fragment_lengths"]["baseline_lengths"] == [10]
    assert body["fragment_lengths"]["variant_lengths"] == [10]
    assert body["fragment_lengths"]["added_lengths"] == []
    assert body["fragment_lengths"]["removed_lengths"] == []


def test_variant_analysis_is_deterministic_and_sorted(client):
    # 两个潜在位点同时被激活，added_cuts 按 (top_cut, bottom_cut) 排序
    sub_id = _create(client, sequence="GAATTTAAAAGAATTT", topology="circular")
    variants = [{"position": 15, "base": "C"}, {"position": 5, "base": "C"}]
    r1 = _analyze(client, sub_id, variants)
    r2 = _analyze(client, sub_id, variants)
    assert r1.status_code == 200 and r2.status_code == 200
    assert r1.json() == r2.json()  # 相同输入，响应完全一致
    keys = [(c["top_cut"], c["bottom_cut"]) for c in r1.json()["added_cuts"]]
    assert keys == [(1, 5), (11, 15)] == sorted(keys)


def test_variant_does_not_modify_saved_record(client):
    sub_id = _create(client)
    before = client.get(f"/digests/{sub_id}").json()
    resp = _analyze(client, sub_id, [{"position": 9, "base": "T"}])
    assert resp.status_code == 200
    assert client.get(f"/digests/{sub_id}").json() == before  # 原记录不变
    assert len(client.get("/digests").json()) == 1  # 未产生新记录


def test_variant_lowercase_base_normalized(client):
    sub_id = _create(client, sequence="AAAAGAATTTTTTT")
    resp = _analyze(client, sub_id, [{"position": 9, "base": "c"}])
    assert resp.status_code == 200
    assert resp.json()["variants_applied"][0]["base"] == "C"


def test_variant_position_out_of_bounds(client):
    sub_id = _create(client)  # 长度 14，合法位置 0..13
    resp = _analyze(client, sub_id, [{"position": 14, "base": "C"}])
    assert resp.status_code == 422
    resp = _analyze(client, sub_id, [{"position": 0, "base": "C"}, {"position": 100, "base": "A"}])
    assert resp.status_code == 422


def test_variant_negative_position_rejected(client):
    sub_id = _create(client)
    assert _analyze(client, sub_id, [{"position": -1, "base": "C"}]).status_code == 422


def test_variant_invalid_base_rejected(client):
    sub_id = _create(client)
    for bad in ("Z", "", "AC", "1"):
        resp = _analyze(client, sub_id, [{"position": 0, "base": bad}])
        assert resp.status_code == 422, bad


def test_variant_duplicate_positions_rejected(client):
    sub_id = _create(client)
    resp = _analyze(client, sub_id, [{"position": 3, "base": "C"},
                                     {"position": 3, "base": "T"}])
    assert resp.status_code == 422


def test_empty_variants_rejected(client):
    sub_id = _create(client)
    assert _analyze(client, sub_id, []).status_code == 422


def test_variant_analysis_missing_record_404(client):
    resp = _analyze(client, 999, [{"position": 0, "base": "C"}])
    assert resp.status_code == 404
