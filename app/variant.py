"""变异影响分析。

对已保存的酶切记录应用一组 0-based 单碱基替换（SNV），在内存中重跑
酶切并与原结果对比，返回新增/消失/末端改变的切口以及片段长度差异。
分析过程只读取记录，不写数据库，原记录保持不变。

切口身份约定
------------
物理切口以 (top_cut, bottom_cut) 为唯一标识。单碱基替换不改变序列长度，
因此同一切口的坐标不会平移：识别位点被破坏/新建只会导致切口消失/新增，
"坐标改变"在结果中表现为一条 removed 加一条 added。若替换发生在黏端
突出序列内部且未破坏识别位点（如模式含简并碱基），同一切口的 overhang
会改变，归入 changed。
"""

from __future__ import annotations

from collections import Counter

from .digest import run_digest


def apply_variants(sequence: str, variants) -> str:
    """按 0-based 位置应用单碱基替换，返回突变后的序列。

    调用方需保证位置合法（0 <= position < len(sequence)）且无重复。
    """
    bases = list(sequence)
    for v in variants:
        bases[v.position] = v.base
    return "".join(bases)


def _cut_key(cut: dict) -> tuple[int, int]:
    return (cut["top_cut"], cut["bottom_cut"])


def diff_cuts(baseline_cuts: list[dict], variant_cuts: list[dict]):
    """对比两次酶切的切口列表，返回 (added, removed, changed)。

    added/removed 为完整切口记录，按 (top_cut, bottom_cut) 排序；
    changed 为同一坐标切口末端特征（end_type/overhang）的变化，
    同样按坐标排序，保证输出顺序稳定。
    """
    before = {_cut_key(c): c for c in baseline_cuts}
    after = {_cut_key(c): c for c in variant_cuts}

    added = sorted((after[k] for k in after.keys() - before.keys()), key=_cut_key)
    removed = sorted((before[k] for k in before.keys() - after.keys()), key=_cut_key)

    changed = []
    for key in sorted(before.keys() & after.keys()):
        old, new = before[key], after[key]
        changes = []
        if old["end_type"] != new["end_type"]:
            changes.append("end_type")
        if old["overhang"] != new["overhang"]:
            changes.append("overhang")
        if changes:
            changed.append(
                {
                    "top_cut": key[0],
                    "bottom_cut": key[1],
                    "changes": changes,
                    "before": {"end_type": old["end_type"], "overhang": old["overhang"]},
                    "after": {"end_type": new["end_type"], "overhang": new["overhang"]},
                }
            )
    return added, removed, changed


def diff_fragment_lengths(baseline_fragments: list[dict], variant_fragments: list[dict]) -> dict:
    """片段长度差异：两边长度序列（按片段顺序）及多重集差（升序）。"""
    before = [f["length"] for f in baseline_fragments]
    after = [f["length"] for f in variant_fragments]
    return {
        "baseline_lengths": before,
        "variant_lengths": after,
        "added_lengths": sorted((Counter(after) - Counter(before)).elements()),
        "removed_lengths": sorted((Counter(before) - Counter(after)).elements()),
    }


def analyze_variants(sequence: str, topology: str, enzymes, variants) -> dict:
    """应用变异并重跑酶切，返回与原始序列结果的完整对比。"""
    mutated = apply_variants(sequence, variants)
    baseline = run_digest(sequence, topology, enzymes)
    variant = run_digest(mutated, topology, enzymes)
    added, removed, changed = diff_cuts(baseline["cuts"], variant["cuts"])
    return {
        "mutated_sequence": mutated,
        "variants_applied": [
            {"position": v.position, "ref_base": sequence[v.position], "base": v.base}
            for v in sorted(variants, key=lambda v: v.position)
        ],
        "added_cuts": added,
        "removed_cuts": removed,
        "changed_cuts": changed,
        "fragment_lengths": diff_fragment_lengths(baseline["fragments"], variant["fragments"]),
    }
