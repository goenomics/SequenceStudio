"""Primer design — Tm, GC%, hairpin, and primer pair generation."""

from __future__ import annotations

from dataclasses import dataclass
from math import log

from Bio.Seq import Seq


@dataclass
class Primer:
    sequence: str
    start: int       # 0-based on template (forward strand)
    end: int         # 0-based exclusive
    strand: int      # +1 forward, -1 reverse
    tm: float
    gc_pct: float
    length: int
    name: str = ""

    @property
    def display_seq(self) -> str:
        """Sequence as it would be ordered in a tube (5'→3')."""
        return self.sequence


# ---------------------------------------------------------------------------
# Thermodynamic helpers
# ---------------------------------------------------------------------------

def gc_content(seq: str) -> float:
    s = seq.upper()
    return ((s.count("G") + s.count("C")) / len(s) * 100) if s else 0.0


# SantaLucia 1998 nearest-neighbor parameters: (ΔH kcal/mol, ΔS cal/mol/K)
_NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "AC": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "TA": (-7.2, -21.3), "TT": (-7.9, -22.2), "TC": (-8.2, -22.2), "TG": (-8.5, -22.7),
    "CA": (-8.5, -22.7), "CT": (-7.8, -21.0), "CC": (-8.0, -19.9), "CG": (-10.6, -27.2),
    "GA": (-8.2, -22.2), "GT": (-7.8, -21.0), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
}
_R = 1.987          # gas constant cal/mol/K
_CT = 250e-9        # typical primer concentration 250 nM


def calculate_tm(seq: str) -> float:
    """Melting temperature using nearest-neighbor thermodynamics (SantaLucia 1998)."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    if n < 14:
        # Wallace rule for very short oligos
        return 2.0 * (s.count("A") + s.count("T")) + 4.0 * (s.count("G") + s.count("C"))

    dH = 0.2   # initiation (kcal/mol)
    dS = -5.7  # initiation (cal/mol/K)
    for i in range(n - 1):
        di = s[i : i + 2]
        if di in _NN:
            h, sv = _NN[di]
            dH += h
            dS += sv

    denom = dS + _R * log(_CT / 4.0)
    if denom == 0:
        return 0.0
    return round((dH * 1000.0) / denom - 273.15, 1)


def has_hairpin(seq: str, min_stem: int = 4, min_loop: int = 3) -> bool:
    """Return True if the sequence can form a simple hairpin."""
    s = seq.upper()
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    n = len(s)
    for stem in range(min_stem, n // 2 + 1):
        for i in range(n - stem * 2 - min_loop + 1):
            s1 = s[i : i + stem]
            s2 = s[i + stem + min_loop : i + stem * 2 + min_loop]
            if len(s2) < stem:
                continue
            rc2 = "".join(comp.get(c, "N") for c in reversed(s2))
            if s1 == rc2:
                return True
    return False


# ---------------------------------------------------------------------------
# Primer design
# ---------------------------------------------------------------------------

def design_primers(
    template: str,
    target_start: int,
    target_end: int,
    length_range: tuple[int, int] = (18, 25),
    tm_range: tuple[float, float] = (55.0, 65.0),
    gc_range: tuple[float, float] = (40.0, 60.0),
    max_per_strand: int = 10,
) -> list[Primer]:
    """
    Design PCR primer pairs flanking *template[target_start:target_end]*.
    Returns the best forward primers first, then best reverse primers.
    """
    tmpl = template.upper().replace("-", "")
    min_len, max_len = length_range
    target_tm = sum(tm_range) / 2.0

    candidates: list[Primer] = []

    # Forward primers ending before target_start
    fwd_search_start = max(0, target_start - 300)
    for pos in range(fwd_search_start, target_start):
        for plen in range(min_len, max_len + 1):
            if pos + plen > target_start:
                break
            seq = tmpl[pos : pos + plen]
            if "N" in seq:
                continue
            tm = calculate_tm(seq)
            gc = gc_content(seq)
            if tm_range[0] <= tm <= tm_range[1] and gc_range[0] <= gc <= gc_range[1]:
                if not has_hairpin(seq):
                    candidates.append(Primer(
                        sequence=seq,
                        start=pos, end=pos + plen,
                        strand=1, tm=tm, gc_pct=round(gc, 1),
                        length=plen, name=f"FWD_{pos + 1}",
                    ))

    # Reverse primers starting after target_end
    rev_search_end = min(len(tmpl), target_end + 300)
    for pos in range(target_end, rev_search_end):
        for plen in range(min_len, max_len + 1):
            if pos + plen > rev_search_end:
                break
            fwd_seq = tmpl[pos : pos + plen]
            if "N" in fwd_seq:
                continue
            rev_seq = str(Seq(fwd_seq).reverse_complement())
            tm = calculate_tm(rev_seq)
            gc = gc_content(rev_seq)
            if tm_range[0] <= tm <= tm_range[1] and gc_range[0] <= gc <= gc_range[1]:
                if not has_hairpin(rev_seq):
                    candidates.append(Primer(
                        sequence=rev_seq,
                        start=pos, end=pos + plen,
                        strand=-1, tm=tm, gc_pct=round(gc, 1),
                        length=plen, name=f"REV_{pos + 1}",
                    ))

    fwd = sorted([p for p in candidates if p.strand == 1],
                 key=lambda p: abs(p.tm - target_tm))
    rev = sorted([p for p in candidates if p.strand == -1],
                 key=lambda p: abs(p.tm - target_tm))

    return fwd[:max_per_strand] + rev[:max_per_strand]
