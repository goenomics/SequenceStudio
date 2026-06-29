"""File I/O — read and write common biological sequence formats."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO, AlignIO

from sequencestudio.models.sequence import SequenceRecord, SequenceAlignment

# Extension → BioPython format name
_READ_FMT: dict[str, str] = {
    ".fasta": "fasta",
    ".fa":    "fasta",
    ".fna":   "fasta",
    ".ffn":   "fasta",
    ".faa":   "fasta",
    ".frn":   "fasta",
    ".gb":    "genbank",
    ".gbk":   "genbank",
    ".genbank": "genbank",
    ".aln":   "clustal",
    ".clustal": "clustal",
    ".nex":   "nexus",
    ".nexus": "nexus",
    ".phy":   "phylip",
    ".phylip": "phylip",
    ".embl":  "embl",
}

WRITE_FORMATS: dict[str, str] = {
    "FASTA":   "fasta",
    "GenBank": "genbank",
    "Clustal": "clustal",
    "NEXUS":   "nexus",
    "PHYLIP":  "phylip",
}


def detect_format(filepath: str) -> str:
    ext = Path(filepath).suffix.lower()
    return _READ_FMT.get(ext, "fasta")


def read_sequences(filepath: str) -> list[SequenceRecord]:
    """Read all sequences from a file, auto-detecting format."""
    fmt = detect_format(filepath)
    records: list[SequenceRecord] = []
    with open(filepath, "r") as fh:
        for rec in SeqIO.parse(fh, fmt):
            records.append(SequenceRecord(record=rec))
    if not records:
        raise IOError(f"No sequences found in {filepath!r} (format detected: {fmt})")
    return records


def write_sequences(filepath: str, records: list[SequenceRecord], fmt: str = "fasta") -> None:
    """Write sequences to a file using the given BioPython format string."""
    bio_records = [r.record for r in records]
    with open(filepath, "w") as fh:
        SeqIO.write(bio_records, fh, fmt)


def read_alignment(filepath: str) -> SequenceAlignment:
    """Read an alignment file and wrap it as a SequenceAlignment."""
    records = read_sequences(filepath)
    return SequenceAlignment(records=records, name=Path(filepath).stem)
