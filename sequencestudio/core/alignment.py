"""Multiple sequence alignment.

External tools (fastest, best quality):
  brew install mafft        ← recommended, always available
  brew install muscle

Built-in fallback:
  Pure-Python progressive pairwise alignment via BioPython PairwiseAligner.
  Suitable for small sets (< ~50 sequences of < ~5 kb each).
"""

from __future__ import annotations

import os
import subprocess
import tempfile

from Bio import SeqIO, AlignIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from sequencestudio.models.sequence import SequenceRecord, SequenceAlignment


def _tool_available(name: str) -> bool:
    try:
        return subprocess.run(["which", name], capture_output=True).returncode == 0
    except Exception:
        return False


# ---------------------------------------------------------------------------
# MAFFT (recommended — brew install mafft)
# ---------------------------------------------------------------------------

def align_with_mafft(records: list[SequenceRecord]) -> SequenceAlignment:
    if not _tool_available("mafft"):
        raise RuntimeError(
            "MAFFT not found.\n"
            "Install with:  brew install mafft"
        )
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "input.fasta")
        out = os.path.join(tmp, "output.fasta")
        with open(inp, "w") as fh:
            for r in records:
                fh.write(f">{r.id}\n{r.sequence}\n")
        proc = subprocess.run(
            ["mafft", "--auto", "--quiet", inp],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"MAFFT error:\n{proc.stderr}")
        # mafft writes to stdout
        with open(out, "w") as fh:
            fh.write(proc.stdout)
        aligned = [SequenceRecord(record=r) for r in SeqIO.parse(out, "fasta")]
    return SequenceAlignment(records=aligned)


# ---------------------------------------------------------------------------
# ClustalOmega (brew install clustal-omega — source sometimes unavailable)
# ---------------------------------------------------------------------------

def align_with_clustalo(records: list[SequenceRecord]) -> SequenceAlignment:
    if not _tool_available("clustalo"):
        raise RuntimeError(
            "ClustalOmega (clustalo) not found.\n"
            "The Homebrew formula source is currently unavailable (HTTP 403).\n"
            "Use MAFFT instead:  brew install mafft"
        )
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "input.fasta")
        out = os.path.join(tmp, "output.aln")
        with open(inp, "w") as fh:
            for r in records:
                fh.write(f">{r.id}\n{r.sequence}\n")
        proc = subprocess.run(
            ["clustalo", "-i", inp, "-o", out, "--outfmt=clustal", "--force"],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"ClustalOmega error:\n{proc.stderr}")
        aligned = [SequenceRecord(record=r) for r in AlignIO.read(out, "clustal")]
    return SequenceAlignment(records=aligned)


# ---------------------------------------------------------------------------
# MUSCLE (brew install muscle)
# ---------------------------------------------------------------------------

def align_with_muscle(records: list[SequenceRecord]) -> SequenceAlignment:
    if not _tool_available("muscle"):
        raise RuntimeError(
            "MUSCLE not found.\n"
            "Install with:  brew install muscle"
        )
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "input.fasta")
        out = os.path.join(tmp, "output.fasta")
        with open(inp, "w") as fh:
            for r in records:
                fh.write(f">{r.id}\n{r.sequence}\n")
        # Try MUSCLE v5 syntax first, fall back to v3
        proc = subprocess.run(
            ["muscle", "-align", inp, "-output", out],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            proc = subprocess.run(
                ["muscle", "-in", inp, "-out", out],
                capture_output=True, text=True,
            )
        if proc.returncode != 0:
            raise RuntimeError(f"MUSCLE error:\n{proc.stderr}")
        aligned = [SequenceRecord(record=r) for r in SeqIO.parse(out, "fasta")]
    return SequenceAlignment(records=aligned)


# ---------------------------------------------------------------------------
# Pure-Python progressive alignment fallback (no external tools required)
# Uses BioPython PairwiseAligner — quality is lower than MAFFT/MUSCLE but
# works for small sets when no external aligner is installed.
# ---------------------------------------------------------------------------

def _pairwise_align(seq_a: str, seq_b: str) -> tuple[str, str]:
    """
    Globally align two sequences and return the gap-padded string pair.
    Uses BioPython PairwiseAligner coordinate blocks to reconstruct gapped seqs.
    """
    clean_a = seq_a.replace("-", "")
    clean_b = seq_b.replace("-", "")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best = next(iter(aligner.align(clean_a, clean_b)))

    # best.aligned → (target_blocks, query_blocks)
    # each block list: array of (start, end) pairs of *ungapped* positions
    t_blocks = best.aligned[0]
    q_blocks  = best.aligned[1]

    result_a: list[str] = []
    result_b: list[str] = []
    prev_t = prev_q = 0

    for (t0, t1), (q0, q1) in zip(t_blocks, q_blocks):
        # gap in a (insert in b before this block)
        if q0 > prev_q:
            result_a.append("-" * (q0 - prev_q))
            result_b.append(clean_b[prev_q:q0])
        # gap in b (insert in a before this block)
        if t0 > prev_t:
            result_a.append(clean_a[prev_t:t0])
            result_b.append("-" * (t0 - prev_t))
        # aligned block
        result_a.append(clean_a[t0:t1])
        result_b.append(clean_b[q0:q1])
        prev_t, prev_q = t1, q1

    # trailing residues
    if prev_t < len(clean_a):
        result_a.append(clean_a[prev_t:])
        result_b.append("-" * (len(clean_a) - prev_t))
    if prev_q < len(clean_b):
        result_b.append(clean_b[prev_q:])
        result_a.append("-" * (len(clean_b) - prev_q))

    a_str = "".join(result_a)
    b_str = "".join(result_b)
    maxl = max(len(a_str), len(b_str))
    return a_str.ljust(maxl, "-"), b_str.ljust(maxl, "-")


def align_builtin(records: list[SequenceRecord]) -> SequenceAlignment:
    """
    Progressive pairwise alignment using BioPython only — no external tool needed.
    Suitable for ≤ 50 sequences of ≤ 5,000 bp each.
    """
    if len(records) < 2:
        raise ValueError("Need at least 2 sequences to align.")

    # Seed with first pair
    a0, a1 = _pairwise_align(records[0].sequence, records[1].sequence)
    aligned: list[str] = [a0, a1]

    # Progressive addition: align each new sequence to the running consensus
    for rec in records[2:]:
        curr_len = len(aligned[0])

        # Plurality-vote consensus (ignoring gaps)
        cons_chars: list[str] = []
        for col in range(curr_len):
            chars = [s[col] for s in aligned]
            non_gap = [c for c in chars if c != "-"]
            cons_chars.append(max(set(non_gap), key=non_gap.count) if non_gap else "-")
        cons_str = "".join(cons_chars)

        cons_aln, new_aln = _pairwise_align(cons_str, rec.sequence)

        # Re-gap existing sequences according to the new gaps introduced in cons_aln
        re_gapped: list[str] = []
        for existing in aligned:
            out: list[str] = []
            ei = 0
            for ch in cons_aln:
                if ch == "-":
                    out.append("-")           # new gap column
                else:
                    out.append(existing[ei] if ei < len(existing) else "-")
                    ei += 1
            re_gapped.append("".join(out))

        maxl = max(len(s) for s in re_gapped + [new_aln])
        aligned = [s.ljust(maxl, "-") for s in re_gapped] + [new_aln.ljust(maxl, "-")]

    # Build SequenceRecord list
    result: list[SequenceRecord] = []
    for rec, aln_seq in zip(records, aligned):
        bio = SeqRecord(Seq(aln_seq), id=rec.id, name=rec.name, description=rec.description)
        result.append(SequenceRecord(record=bio))

    return SequenceAlignment(records=result)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def align_sequences(records: list[SequenceRecord], method: str = "mafft") -> SequenceAlignment:
    """
    Align *records* using the specified *method*.

    method values:
      'mafft'    — MAFFT (brew install mafft)          ← recommended
      'muscle'   — MUSCLE (brew install muscle)
      'clustalo' — ClustalOmega (source currently 403; use mafft instead)
      'builtin'  — Pure-Python progressive fallback (no install required)
    """
    if method == "mafft":
        return align_with_mafft(records)
    elif method == "muscle":
        return align_with_muscle(records)
    elif method == "clustalo":
        return align_with_clustalo(records)
    elif method == "builtin":
        return align_builtin(records)
    raise ValueError(f"Unknown alignment method: {method!r}")
