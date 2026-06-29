"""Core sequence and alignment data models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from Bio.SeqRecord import SeqRecord


@dataclass
class SequenceRecord:
    """UI wrapper around a BioPython SeqRecord, carrying display metadata."""

    record: SeqRecord
    selected: bool = False
    group: Optional[str] = None

    # ------------------------------------------------------------------
    # Properties delegating to the underlying SeqRecord
    # ------------------------------------------------------------------

    @property
    def id(self) -> str:
        return self.record.id

    @property
    def name(self) -> str:
        return self.record.name

    @property
    def description(self) -> str:
        return self.record.description

    @property
    def sequence(self) -> str:
        return str(self.record.seq)

    @property
    def length(self) -> int:
        return len(self.record.seq)

    @property
    def molecule_type(self) -> str:
        """Heuristic detection of DNA / RNA / Protein."""
        raw = self.sequence.upper().replace("-", "").replace(".", "").replace("N", "")
        if not raw:
            return "Unknown"
        chars = set(raw)
        if chars.issubset(set("ACGT") | set("RYSWKMBDHV")):
            return "DNA"
        if chars.issubset(set("ACGU") | set("N")):
            return "RNA"
        if chars.issubset(set("ACDEFGHIKLMNPQRSTVWY") | set("*XBZ")):
            return "Protein"
        return "DNA"  # conservative default

    def gc_content(self) -> float:
        """GC percentage, ignoring gaps."""
        seq = self.sequence.upper()
        gc = seq.count("G") + seq.count("C")
        total = len(seq.replace("-", "").replace(".", ""))
        return (gc / total * 100) if total else 0.0


@dataclass
class SequenceAlignment:
    """A multiple sequence alignment — list of same-length SequenceRecords."""

    records: list[SequenceRecord] = field(default_factory=list)
    name: str = "Untitled"

    @property
    def num_sequences(self) -> int:
        return len(self.records)

    @property
    def alignment_length(self) -> int:
        if not self.records:
            return 0
        return max(r.length for r in self.records)

    def add_record(self, record: SequenceRecord) -> None:
        self.records.append(record)

    def remove_record(self, index: int) -> None:
        if 0 <= index < len(self.records):
            self.records.pop(index)

    def consensus(self) -> str:
        """Plurality-vote consensus across alignment columns."""
        if not self.records:
            return ""
        result: list[str] = []
        for col in range(self.alignment_length):
            counts: dict[str, int] = {}
            for rec in self.records:
                ch = rec.sequence[col].upper() if col < rec.length else "-"
                if ch not in ("-", "."):
                    counts[ch] = counts.get(ch, 0) + 1
            result.append(max(counts, key=counts.get) if counts else "-")
        return "".join(result)

    def conservation(self) -> list[float]:
        """Per-column conservation score in [0, 1]."""
        if not self.records:
            return []
        n = self.num_sequences
        result: list[float] = []
        for col in range(self.alignment_length):
            counts: dict[str, int] = {}
            for rec in self.records:
                ch = rec.sequence[col].upper() if col < rec.length else "-"
                counts[ch] = counts.get(ch, 0) + 1
            max_count = max(counts.values()) if counts else 0
            result.append(max_count / n if n else 0.0)
        return result
