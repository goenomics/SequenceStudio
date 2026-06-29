"""Open Reading Frame (ORF) finder — all 6 reading frames."""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq


@dataclass
class ORF:
    frame: int      # +1/+2/+3 for forward, -1/-2/-3 for reverse
    start: int      # 0-based on the *original* (forward) sequence
    end: int        # 0-based exclusive on the original sequence
    length: int     # nucleotide length
    protein: str    # amino-acid translation (stop codon excluded)


def find_orfs(sequence: str, min_nt_length: int = 150) -> list[ORF]:
    """
    Scan all 6 reading frames of *sequence* for ORFs ≥ *min_nt_length* nt.
    Returns a list sorted by length (longest first).
    """
    seq = sequence.upper().replace("-", "").replace(".", "")
    fwd = Seq(seq)
    rev = fwd.reverse_complement()
    seq_len = len(seq)
    results: list[ORF] = []

    def _scan(dna: Seq, strand: int, frame_offset: int) -> None:
        aa_seq = str(dna[frame_offset:].translate(to_stop=False))
        in_orf = False
        orf_aa_start = 0

        for i, aa in enumerate(aa_seq):
            if aa == "M" and not in_orf:
                in_orf = True
                orf_aa_start = i
            if in_orf and aa == "*":
                nt_len = (i - orf_aa_start) * 3
                if nt_len >= min_nt_length:
                    protein = aa_seq[orf_aa_start:i]
                    nt_start_in_frame = frame_offset + orf_aa_start * 3
                    nt_end_in_frame   = frame_offset + i * 3

                    if strand == 1:
                        results.append(ORF(
                            frame=frame_offset + 1,
                            start=nt_start_in_frame,
                            end=nt_end_in_frame,
                            length=nt_len,
                            protein=protein,
                        ))
                    else:
                        # Map reverse-complement coordinates back to forward strand
                        fwd_start = seq_len - nt_end_in_frame
                        fwd_end   = seq_len - nt_start_in_frame
                        results.append(ORF(
                            frame=-(frame_offset + 1),
                            start=fwd_start,
                            end=fwd_end,
                            length=nt_len,
                            protein=protein,
                        ))
                in_orf = False

    for offset in range(3):
        _scan(fwd, 1, offset)
        _scan(rev, -1, offset)

    return sorted(results, key=lambda o: o.length, reverse=True)
