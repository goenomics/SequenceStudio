"""Gene prediction backends (pyrodigal + ORF fallback)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from Bio.Seq import Seq

from sequencestudio.core.orf import find_orfs


@dataclass
class GenePrediction:
    """Normalized gene feature output."""

    start: int            # 0-based, inclusive
    end: int              # 0-based, exclusive
    strand: str           # '+' or '-'
    source: str           # predictor name
    length_nt: int
    feature_type: str = "CDS"
    feature_id: str = ""
    parent_id: str = ""
    phase: Optional[int] = None
    protein: str = ""
    confidence: Optional[float] = None


@dataclass
class PredictionRun:
    """Result payload including metadata about the selected backend."""

    method: str
    notes: str
    genes: list[GenePrediction]


@dataclass
class ExportSequence:
    """Derived transcript/protein sequence for FASTA export."""

    record_id: str
    sequence: str
    description: str = ""


def import_predictions_gff3(text: str, preferred_seq_id: Optional[str] = None) -> tuple[str, PredictionRun]:
    """Parse a GFF3 document into typed feature records.

    Supports `gene`, `mRNA`, `exon`, and `CDS` features and preserves IDs/parents.
    Returns `(seq_id, PredictionRun)`.
    """
    seq_features: dict[str, list[GenePrediction]] = {}

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")
        if len(parts) != 9:
            continue

        seq_id, source, feature_type, start_s, end_s, score_s, strand, _phase, attrs = parts
        if feature_type not in {"gene", "mRNA", "exon", "CDS"}:
            continue

        try:
            start0 = max(0, int(start_s) - 1)
            end0 = max(start0, int(end_s))
        except ValueError:
            continue

        attr_map: dict[str, str] = {}
        for item in attrs.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                attr_map[k] = v

        confidence = None if score_s == "." else float(score_s)
        phase = None if _phase == "." else int(_phase)
        pred = GenePrediction(
            start=start0,
            end=end0,
            strand=strand if strand in {"+", "-"} else "+",
            source=source or "gff3",
            length_nt=end0 - start0,
            feature_type=feature_type,
            feature_id=attr_map.get("ID", ""),
            parent_id=attr_map.get("Parent", ""),
            phase=phase,
            protein="",
            confidence=confidence,
        )
        seq_features.setdefault(seq_id, []).append(pred)

    available_ids = set(seq_features)
    if not available_ids:
        raise ValueError("No supported gene features found in GFF3 input")

    if preferred_seq_id and preferred_seq_id in available_ids:
        seq_id = preferred_seq_id
    else:
        seq_id = next(iter(available_ids))

    genes = sorted(
        seq_features.get(seq_id, []),
        key=lambda g: (g.start, g.end, g.strand, g.feature_type),
    )
    return seq_id, PredictionRun(
        method="gff3-import",
        notes="Imported from GFF3 annotation file with typed features.",
        genes=genes,
    )


def export_predictions_gff3(run: PredictionRun, seq_id: str) -> str:
    """Serialize predictions to a minimal GFF3 document."""
    lines = ["##gff-version 3"]
    for idx, gene in enumerate(run.genes, start=1):
        gene_id = gene.feature_id or f"gene{idx}"
        cds_id = f"cds{idx}"
        score = "." if gene.confidence is None else f"{gene.confidence:.3f}"
        start1 = gene.start + 1
        end1 = gene.end
        attrs_gene = f"ID={gene_id};Name={gene_id};source={gene.source}"
        attrs_cds = f"ID={cds_id};Parent={gene_id};product=predicted_protein"

        if gene.feature_type in {"gene", "mRNA", "exon", "CDS"} and gene.feature_id:
            attrs = [f"ID={gene.feature_id}"]
            if gene.parent_id:
                attrs.append(f"Parent={gene.parent_id}")
            attrs.append(f"source={gene.source}")
            lines.append(
                "\t".join([
                    seq_id,
                    run.method,
                    gene.feature_type,
                    str(start1),
                    str(end1),
                    score,
                    gene.strand,
                    "0" if gene.feature_type == "CDS" else ".",
                    ";".join(attrs),
                ])
            )
            continue

        lines.append(
            "\t".join([
                seq_id,
                run.method,
                "gene",
                str(start1),
                str(end1),
                score,
                gene.strand,
                ".",
                attrs_gene,
            ])
        )
        lines.append(
            "\t".join([
                seq_id,
                run.method,
                "CDS",
                str(start1),
                str(end1),
                score,
                gene.strand,
                "0",
                attrs_cds,
            ])
        )
    return "\n".join(lines) + "\n"


def _clean_sequence(sequence: str) -> str:
    return sequence.upper().replace("-", "").replace(".", "")


def _predict_with_orf_fallback(sequence: str, min_nt_length: int) -> list[GenePrediction]:
    genes: list[GenePrediction] = []
    for orf in find_orfs(sequence, min_nt_length=min_nt_length):
        genes.append(
            GenePrediction(
                start=orf.start,
                end=orf.end,
                strand="+" if orf.frame > 0 else "-",
                source="orf-fallback",
                length_nt=orf.length,
                feature_type="CDS",
                feature_id="",
                parent_id="",
                phase=0,
                protein=orf.protein,
                confidence=None,
            )
        )
    genes.sort(key=lambda g: (g.start, g.end))
    return genes


def _predict_with_pyrodigal(sequence: str, min_nt_length: int) -> list[GenePrediction]:
    import pyrodigal  # optional dependency

    dna = _clean_sequence(sequence)
    if not dna:
        return []

    finder = pyrodigal.GeneFinder(meta=True)
    # pyrodigal accepts DNA sequence strings and returns predictions with 1-based coordinates
    results = finder.find_genes(dna)

    genes: list[GenePrediction] = []
    for pred in results:
        begin = int(getattr(pred, "begin", 1))
        end = int(getattr(pred, "end", begin))
        strand_val = int(getattr(pred, "strand", 1))

        start0 = max(0, begin - 1)
        end0 = max(start0, end)
        length_nt = end0 - start0
        if length_nt < min_nt_length:
            continue

        prot = ""
        translate_fn = getattr(pred, "translate", None)
        if callable(translate_fn):
            try:
                prot = str(translate_fn())
            except Exception:
                prot = ""

        score_val = getattr(pred, "score", None)
        confidence = float(score_val) if score_val is not None else None

        genes.append(
            GenePrediction(
                start=start0,
                end=end0,
                strand="+" if strand_val >= 0 else "-",
                source="pyrodigal",
                length_nt=length_nt,
                feature_type="CDS",
                feature_id="",
                parent_id="",
                phase=0,
                protein=prot,
                confidence=confidence,
            )
        )

    genes.sort(key=lambda g: (g.start, g.end))
    return genes


def predict_genes(sequence: str, method: str = "auto", min_nt_length: int = 150) -> PredictionRun:
    """Predict genes from DNA sequence.

        Methods:
        - auto: prefer pyrodigal, otherwise ORF fallback
        - pyrodigal: force pyrodigal backend
        - orf: force built-in ORF fallback
        - eukaryote: reserved backend slot for future AUGUSTUS/BRAKER integration;
            currently falls back to ORF scan with an explanatory note
    """
    dna = _clean_sequence(sequence)
    if not dna:
        return PredictionRun(method="none", notes="Empty sequence.", genes=[])

    if method not in {"auto", "pyrodigal", "orf", "eukaryote"}:
        raise ValueError(f"Unsupported method: {method}")

    if method == "eukaryote":
        genes = _predict_with_orf_fallback(dna, min_nt_length=min_nt_length)
        return PredictionRun(
            method="eukaryote-slot",
            notes=(
                "Future backend slot for AUGUSTUS/BRAKER-style eukaryotic prediction; "
                "currently using ORF fallback."
            ),
            genes=genes,
        )

    if method == "orf":
        genes = _predict_with_orf_fallback(dna, min_nt_length=min_nt_length)
        return PredictionRun(
            method="orf",
            notes="Built-in ORF fallback (6-frame scan).",
            genes=genes,
        )

    try:
        genes = _predict_with_pyrodigal(dna, min_nt_length=min_nt_length)
        return PredictionRun(
            method="pyrodigal",
            notes="Predicted with pyrodigal (prokaryotic gene finder).",
            genes=genes,
        )
    except Exception as exc:
        if method == "pyrodigal":
            raise RuntimeError(f"pyrodigal backend failed: {exc}") from exc

        genes = _predict_with_orf_fallback(dna, min_nt_length=min_nt_length)
        return PredictionRun(
            method="orf",
            notes=f"pyrodigal unavailable ({exc.__class__.__name__}); using ORF fallback.",
            genes=genes,
        )


def _extract_segments(genome_sequence: str, features: list[GenePrediction]) -> str:
    dna = _clean_sequence(genome_sequence)
    if not features:
        return ""

    strand = features[0].strand
    ordered = sorted(features, key=lambda g: g.start, reverse=(strand == "-"))
    parts = [dna[f.start:f.end] for f in ordered]
    joined = "".join(parts)
    if strand == "-":
        joined = str(Seq(joined).reverse_complement())
    return joined


def derive_transcript_sequences(run: PredictionRun, genome_sequence: str) -> list[ExportSequence]:
    """Derive transcript-like nucleotide sequences from imported/predicted features."""
    genes = run.genes
    exon_groups: dict[str, list[GenePrediction]] = {}
    mrna_features: dict[str, GenePrediction] = {}

    for feature in genes:
        if feature.feature_type == "mRNA" and feature.feature_id:
            mrna_features[feature.feature_id] = feature
        elif feature.feature_type == "exon" and feature.parent_id:
            exon_groups.setdefault(feature.parent_id, []).append(feature)

    exports: list[ExportSequence] = []
    if exon_groups:
        for parent_id, exons in sorted(exon_groups.items()):
            seq = _extract_segments(genome_sequence, exons)
            if not seq:
                continue
            parent = mrna_features.get(parent_id)
            desc = f"type=transcript parent={parent.parent_id if parent else ''}".strip()
            exports.append(ExportSequence(record_id=parent_id, sequence=seq, description=desc))
        return exports

    # Fallback: export mRNA spans, then gene spans, then CDS spans
    fallback_types = ["mRNA", "gene", "CDS"]
    for feature_type in fallback_types:
        subset = [g for g in genes if g.feature_type == feature_type]
        if subset:
            for idx, feature in enumerate(subset, start=1):
                seq = _extract_segments(genome_sequence, [feature])
                record_id = feature.feature_id or f"{feature_type.lower()}{idx}"
                exports.append(
                    ExportSequence(
                        record_id=record_id,
                        sequence=seq,
                        description=f"type={feature_type} strand={feature.strand}",
                    )
                )
            return exports

    return exports


def derive_protein_sequences(run: PredictionRun, genome_sequence: str) -> list[ExportSequence]:
    """Derive protein FASTA records from CDS features or built-in protein strings."""
    genes = run.genes
    cds_groups: dict[str, list[GenePrediction]] = {}
    standalone_cdss: list[GenePrediction] = []

    for idx, feature in enumerate(genes, start=1):
        if feature.feature_type != "CDS":
            continue
        group_id = feature.parent_id or feature.feature_id
        if group_id:
            cds_groups.setdefault(group_id, []).append(feature)
        else:
            standalone_cdss.append(feature)

    exports: list[ExportSequence] = []
    for group_id, cds_list in sorted(cds_groups.items()):
        seq = _extract_segments(genome_sequence, cds_list)
        if not seq:
            continue
        trim = len(seq) % 3
        if trim:
            seq = seq[:-trim]
        protein = str(Seq(seq).translate(to_stop=False)) if seq else ""
        if protein.endswith("*"):
            protein = protein[:-1]
        exports.append(
            ExportSequence(
                record_id=group_id,
                sequence=protein,
                description="type=protein source=CDS",
            )
        )

    for idx, feature in enumerate(standalone_cdss, start=1):
        if feature.protein:
            protein = feature.protein.rstrip("*")
        else:
            seq = _extract_segments(genome_sequence, [feature])
            trim = len(seq) % 3
            if trim:
                seq = seq[:-trim]
            protein = str(Seq(seq).translate(to_stop=False)) if seq else ""
            if protein.endswith("*"):
                protein = protein[:-1]
        exports.append(
            ExportSequence(
                record_id=feature.feature_id or f"cds{idx}",
                sequence=protein,
                description=f"type=protein strand={feature.strand}",
            )
        )

    return exports


def export_fasta(records: list[ExportSequence]) -> str:
    """Serialize export records to FASTA text."""
    lines: list[str] = []
    for rec in records:
        header = rec.record_id if not rec.description else f"{rec.record_id} {rec.description}"
        lines.append(f">{header}")
        for i in range(0, len(rec.sequence), 70):
            lines.append(rec.sequence[i:i + 70])
    return "\n".join(lines) + ("\n" if lines else "")
