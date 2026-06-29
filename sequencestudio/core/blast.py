"""NCBI remote BLAST integration via BioPython."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

from Bio.Blast import NCBIWWW, NCBIXML

BLAST_PROGRAMS = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]

BLAST_DATABASES: dict[str, list[str]] = {
    "blastn":  ["nt", "refseq_rna", "refseq_genomic", "core_nt"],
    "blastp":  ["nr", "refseq_protein", "swissprot", "pdb"],
    "blastx":  ["nr", "refseq_protein", "swissprot"],
    "tblastn": ["nt", "refseq_rna"],
    "tblastx": ["nt", "refseq_rna"],
}


@dataclass
class BlastHit:
    title: str
    accession: str
    length: int
    score: float
    evalue: float
    identity: int
    identity_pct: float
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    alignment_text: str = ""


@dataclass
class BlastResult:
    query_id: str
    program: str
    database: str
    hits: list[BlastHit] = field(default_factory=list)

    @property
    def num_hits(self) -> int:
        return len(self.hits)


def run_blast_ncbi(
    sequence: str,
    seq_id: str,
    program: str = "blastn",
    database: str = "nt",
    hitlist_size: int = 50,
    progress: Optional[Callable[[str], None]] = None,
) -> BlastResult:
    """Submit a remote NCBI BLAST query and return parsed results."""
    if progress:
        progress(f"Submitting {program} against {database} …")

    handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=sequence,
        hitlist_size=hitlist_size,
    )

    if progress:
        progress("Parsing BLAST results …")

    blast_record = next(NCBIXML.parse(handle))

    hits: list[BlastHit] = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            pct = (hsp.identities / hsp.align_length * 100) if hsp.align_length else 0.0
            hits.append(BlastHit(
                title=alignment.title[:120],
                accession=alignment.accession,
                length=alignment.length,
                score=hsp.score,
                evalue=hsp.expect,
                identity=hsp.identities,
                identity_pct=round(pct, 1),
                query_start=hsp.query_start,
                query_end=hsp.query_end,
                subject_start=hsp.sbjct_start,
                subject_end=hsp.sbjct_end,
                alignment_text=(
                    f"Query:   {hsp.query}\n"
                    f"         {hsp.match}\n"
                    f"Subject: {hsp.sbjct}"
                ),
            ))

    return BlastResult(query_id=seq_id, program=program, database=database, hits=hits)
