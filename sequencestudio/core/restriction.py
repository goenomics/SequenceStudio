"""Restriction enzyme analysis using BioPython."""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Restriction import CommOnly, Analysis, RestrictionBatch
from Bio.Seq import Seq


@dataclass
class RestrictionResult:
    enzyme: str
    cut_positions: list[int]          # 1-based cut positions
    num_cuts: int
    fragment_sizes: list[int]         # sorted descending


def get_common_enzymes() -> list[str]:
    """Names of commercially available restriction enzymes."""
    return sorted(str(e) for e in CommOnly)


def analyze_sequence(
    sequence: str,
    enzymes: list[str] | None = None,
    linear: bool = True,
) -> list[RestrictionResult]:
    """
    Find restriction sites for *enzymes* in *sequence*.
    Returns results sorted by cut count (fewest first).
    If *enzymes* is None, the full commercial set is used.
    """
    clean = sequence.upper().replace("-", "").replace(".", "")
    seq = Seq(clean)

    if enzymes:
        try:
            rb = RestrictionBatch(enzymes)
        except Exception:
            rb = CommOnly
    else:
        rb = CommOnly

    analysis = Analysis(rb, seq, linear=linear)
    raw: dict = analysis.full()

    results: list[RestrictionResult] = []
    seq_len = len(clean)

    for enzyme, positions in raw.items():
        if not positions:
            continue
        cuts = sorted(int(p) for p in positions)

        # Fragment size calculation
        if linear:
            boundaries = [0] + cuts + [seq_len]
        else:
            boundaries = cuts if cuts else [0]
        fragments = sorted(
            [boundaries[i + 1] - boundaries[i] for i in range(len(boundaries) - 1)],
            reverse=True,
        )

        results.append(RestrictionResult(
            enzyme=str(enzyme),
            cut_positions=cuts,
            num_cuts=len(cuts),
            fragment_sizes=fragments,
        ))

    return sorted(results, key=lambda r: r.num_cuts)
