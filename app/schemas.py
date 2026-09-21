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


class VariantIn(BaseModel):
    """一个 0-based 单碱基替换：将 position 处的碱基替换为 base。"""

    position: int = Field(ge=0, description="0-based 替换位置", examples=[9])
    base: str = Field(description="替换后的碱基（IUPAC 单字符，大小写均可）", examples=["C"])

    @field_validator("base")
    @classmethod
    def _base_ok(cls, v: str) -> str:
        v = v.strip().upper()
        if len(v) != 1 or v not in VALID_BASES:
            raise ValueError(
                f"variant base must be a single IUPAC base, got {v!r}; "
                "allowed: ACGTRYSWKMBDHVN"
            )
        return v


class VariantAnalysisRequest(BaseModel):
    """变异影响分析请求：一组 0-based 单碱基替换，位置不得重复。"""

    variants: list[VariantIn] = Field(min_length=1)

    @model_validator(mode="after")
    def _positions_unique(self) -> "VariantAnalysisRequest":
        seen: set[int] = set()
        dups: set[int] = set()
        for v in self.variants:
            if v.position in seen:
                dups.add(v.position)
            seen.add(v.position)
        if dups:
            raise ValueError(
                "duplicate variant positions: "
                + ", ".join(str(p) for p in sorted(dups))
            )
        return self


class VariantApplied(BaseModel):
    position: int
    ref_base: str
    base: str


class CutEndState(BaseModel):
    end_type: Literal["blunt", "5prime", "3prime"]
    overhang: str


class ChangedCutOut(BaseModel):
    """同一坐标切口在变异前后的末端特征变化。"""

    top_cut: int
    bottom_cut: int
    changes: list[str]
    before: CutEndState
    after: CutEndState


class FragmentLengthDiff(BaseModel):
    baseline_lengths: list[int]
    variant_lengths: list[int]
    added_lengths: list[int]
    removed_lengths: list[int]


class VariantAnalysisResponse(BaseModel):
    submission_id: int
    topology: Literal["linear", "circular"]
    sequence_length: int
    variants_applied: list[VariantApplied]
    mutated_sequence: str
    added_cuts: list[CutOut]
    removed_cuts: list[CutOut]
    changed_cuts: list[ChangedCutOut]
    fragment_lengths: FragmentLengthDiff
