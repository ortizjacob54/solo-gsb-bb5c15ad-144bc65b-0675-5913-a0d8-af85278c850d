"""请求/响应模型与输入校验。

所有校验在数据库写入之前完成：任何一个酶的偏移非法都会使
整个请求以 422 被拒绝，不会有任何数据落库。
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from .iupac import VALID_BASES


def _check_iupac(value: str, what: str) -> str:
    seq = value.strip().upper()
    if not seq:
        raise ValueError(f"{what} must not be empty")
    bad = sorted(set(seq) - VALID_BASES)
    if bad:
        raise ValueError(
            f"{what} contains non-IUPAC bases: {''.join(bad)}; "
            "allowed: ACGTRYSWKMBDHVN"
        )
    return seq


class EnzymeIn(BaseModel):
    """一种限制性内切酶：识别模式 + 正/反链切口偏移。

    偏移以识别位点 5' 端（top 链坐标）为基准，合法范围为
    [0, len(pattern)]；超出范围（含负数）即非法，整个请求被拒绝。
    """

    name: str = Field(min_length=1, max_length=100, examples=["EcoRI"])
    pattern: str = Field(examples=["GAATTC"])
    cut_top: int = Field(description="top 链切口偏移，0..len(pattern)", examples=[1])
    cut_bottom: int = Field(description="bottom 链切口偏移，0..len(pattern)", examples=[5])

    @field_validator("name")
    @classmethod
    def _name_ok(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("enzyme name must not be blank")
        return v

    @field_validator("pattern")
    @classmethod
    def _pattern_ok(cls, v: str) -> str:
        return _check_iupac(v, "pattern")

    @model_validator(mode="after")
    def _offsets_ok(self) -> "EnzymeIn":
        k = len(self.pattern)
        for field_name in ("cut_top", "cut_bottom"):
            value = getattr(self, field_name)
            if not 0 <= value <= k:
                raise ValueError(
                    f"enzyme '{self.name}': {field_name}={value} is out of "
                    f"range [0, {k}] for pattern '{self.pattern}'"
                )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(examples=["AAAAGAATTCTTTT"])
    topology: Literal["linear", "circular"] = "linear"
    enzymes: list[EnzymeIn] = Field(min_length=1)

    @field_validator("sequence")
    @classmethod
    def _sequence_ok(cls, v: str) -> str:
        return _check_iupac(v, "sequence")


class SiteOut(BaseModel):
    enzyme: str
    start: int
    strand: Literal["+", "-"]
    matched: str


class CutOut(BaseModel):
    top_cut: int
    bottom_cut: int
    end_type: Literal["blunt", "5prime", "3prime"]
    overhang: str
    enzymes: list[str]
    sites: list[SiteOut]


class EndOut(BaseModel):
    type: Literal["blunt", "5prime", "3prime", "none"]
    overhang: str


class FragmentOut(BaseModel):
    index: int
    start: int
    end: int
    length: int
    sequence: str
    wraps_origin: bool
    left_end: EndOut
    right_end: EndOut


class DigestResult(BaseModel):
    cuts: list[CutOut]
    fragments: list[FragmentOut]


class DigestResponse(DigestResult):
    id: int
    topology: Literal["linear", "circular"]
    sequence: str
    sequence_length: int
    enzymes: list[EnzymeIn]
    created_at: str


class DigestSummary(BaseModel):
    id: int
    topology: str
    sequence_length: int
    enzyme_count: int
    cut_count: int
    fragment_count: int
    created_at: str
