"""Gene prediction view with export and track visualization."""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QCheckBox,
    QComboBox,
    QSpinBox,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QFileDialog,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as mpatches

from sequencestudio.models.sequence import SequenceRecord
from sequencestudio.core.gene_prediction import (
    predict_genes,
    export_predictions_gff3,
    import_predictions_gff3,
    derive_transcript_sequences,
    derive_protein_sequences,
    export_fasta,
    PredictionRun,
    GenePrediction,
)

_TRACK_COLOR = {
    "+": "#34AADC",
    "-": "#FF6B6B",
}

_FEATURE_Y = {
    ("gene", "+"): 3.1,
    ("mRNA", "+"): 2.3,
    ("exon", "+"): 1.45,
    ("CDS", "+"): 0.6,
    ("gene", "-"): -0.2,
    ("mRNA", "-"): -1.0,
    ("exon", "-"): -1.85,
    ("CDS", "-"): -2.7,
}

_FEATURE_HEIGHT = {
    "gene": 0.22,
    "mRNA": 0.22,
    "exon": 0.34,
    "CDS": 0.48,
}


class GenePredictionView(QWidget):
    """Simple table/detail UI for gene predictions."""

    def __init__(self, record: SequenceRecord, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._record = record
        self._run_result: Optional[PredictionRun] = None
        self._selected_index: Optional[int] = None
        self._visible_indices: list[int] = []
        self._build_ui()
        self._run_prediction()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel(f"<b>{self._record.id}</b> ({self._record.length:,} bp)"))
        ctrl.addStretch()
        ctrl.addWidget(QLabel("Method:"))

        self._method_cb = QComboBox()
        self._method_cb.addItem("Auto (pyrodigal -> ORF fallback)", "auto")
        self._method_cb.addItem("pyrodigal", "pyrodigal")
        self._method_cb.addItem("Built-in ORF fallback", "orf")
        self._method_cb.addItem("Eukaryote backend slot (future)", "eukaryote")
        ctrl.addWidget(self._method_cb)

        ctrl.addWidget(QLabel("Min length (nt):"))
        self._min_spin = QSpinBox()
        self._min_spin.setRange(60, 10000)
        self._min_spin.setSingleStep(30)
        self._min_spin.setValue(150)
        ctrl.addWidget(self._min_spin)

        run_btn = QPushButton("Predict Genes")
        run_btn.clicked.connect(self._run_prediction)
        ctrl.addWidget(run_btn)

        export_btn = QPushButton("Export GFF3")
        export_btn.clicked.connect(self._export_gff3)
        ctrl.addWidget(export_btn)

        import_btn = QPushButton("Load GFF3")
        import_btn.clicked.connect(self._import_gff3)
        ctrl.addWidget(import_btn)

        export_tx_btn = QPushButton("Export Transcripts")
        export_tx_btn.clicked.connect(self._export_transcripts)
        ctrl.addWidget(export_tx_btn)

        export_prot_btn = QPushButton("Export Proteins")
        export_prot_btn.clicked.connect(self._export_proteins)
        ctrl.addWidget(export_prot_btn)

        root.addLayout(ctrl)

        filters = QHBoxLayout()
        filters.addWidget(QLabel("Show:"))
        self._show_gene = QCheckBox("gene")
        self._show_mrna = QCheckBox("mRNA")
        self._show_exon = QCheckBox("exon")
        self._show_cds = QCheckBox("CDS")
        for widget in (self._show_gene, self._show_mrna, self._show_exon, self._show_cds):
            widget.setChecked(True)
            widget.toggled.connect(self._refresh_visible_features)
            filters.addWidget(widget)
        filters.addStretch()
        root.addLayout(filters)

        self._status = QLabel("Ready")
        self._status.setObjectName("status_label")
        root.addWidget(self._status)

        split = QSplitter(Qt.Orientation.Vertical)

        track_w = QWidget()
        track_l = QVBoxLayout(track_w)
        track_l.setContentsMargins(0, 0, 0, 0)
        self._fig = Figure(figsize=(12, 2.2), tight_layout=True)
        self._canvas = FigureCanvas(self._fig)
        track_l.addWidget(self._canvas)
        split.addWidget(track_w)

        self._table = QTableWidget(0, 8)
        self._table.setHorizontalHeaderLabels([
            "Type", "Start", "End", "Length", "Strand", "Source", "ID", "Parent"
        ])
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.itemSelectionChanged.connect(self._on_row_selected)
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setStretchLastSection(True)
        split.addWidget(self._table)

        self._details = QTextEdit()
        self._details.setReadOnly(True)
        self._details.setPlaceholderText("Select a predicted gene to inspect details.")
        split.addWidget(self._details)
        split.setSizes([180, 320, 180])

        root.addWidget(split)

    def _run_prediction(self) -> None:
        method = str(self._method_cb.currentData())
        min_nt = self._min_spin.value()
        try:
            self._run_result = predict_genes(self._record.sequence, method=method, min_nt_length=min_nt)
        except Exception as exc:
            self._run_result = None
            self._selected_index = None
            self._table.setRowCount(0)
            self._status.setText(f"Prediction failed: {exc}")
            self._details.clear()
            self._draw_track()
            return

        genes = self._run_result.genes
        self._selected_index = 0 if genes else None
        self._status.setText(
            f"{len(genes)} prediction(s)  |  backend: {self._run_result.method}  |  {self._run_result.notes}"
        )

        self._refresh_visible_features()

    def _feature_visible(self, gene: GenePrediction) -> bool:
        return (
            (gene.feature_type == "gene" and self._show_gene.isChecked())
            or (gene.feature_type == "mRNA" and self._show_mrna.isChecked())
            or (gene.feature_type == "exon" and self._show_exon.isChecked())
            or (gene.feature_type == "CDS" and self._show_cds.isChecked())
        )

    def _refresh_visible_features(self) -> None:
        if not self._run_result:
            self._visible_indices = []
            self._populate_table([])
            return

        self._visible_indices = [
            idx for idx, gene in enumerate(self._run_result.genes)
            if self._feature_visible(gene)
        ]

        if self._selected_index not in self._visible_indices:
            self._selected_index = self._visible_indices[0] if self._visible_indices else None

        genes = [self._run_result.genes[idx] for idx in self._visible_indices]
        self._populate_table(genes)

        if self._selected_index in self._visible_indices:
            row = self._visible_indices.index(self._selected_index)
            self._table.selectRow(row)
        elif not genes:
            self._details.clear()

    def _populate_table(self, genes: list[GenePrediction]) -> None:
        self._table.clearContents()
        self._table.setRowCount(len(genes))

        for row, gene in enumerate(genes):
            self._table.setItem(row, 0, QTableWidgetItem(gene.feature_type))
            self._table.setItem(row, 1, QTableWidgetItem(str(gene.start + 1)))
            self._table.setItem(row, 2, QTableWidgetItem(str(gene.end)))
            self._table.setItem(row, 3, QTableWidgetItem(str(gene.length_nt)))
            self._table.setItem(row, 4, QTableWidgetItem(gene.strand))
            self._table.setItem(row, 5, QTableWidgetItem(gene.source))
            self._table.setItem(row, 6, QTableWidgetItem(gene.feature_id))
            self._table.setItem(row, 7, QTableWidgetItem(gene.parent_id))

        if genes:
            self._table.selectRow(0)
        else:
            self._details.clear()
        self._draw_track()

    def _export_transcripts(self) -> None:
        if not self._run_result:
            return
        records = derive_transcript_sequences(self._run_result, self._record.sequence)
        if not records:
            self._status.setText("No transcript-like features available for export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Transcript FASTA",
            f"{self._record.id}.transcripts.fasta",
            "FASTA (*.fasta *.fa *.fna);;All Files (*)",
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(export_fasta(records))
        self._status.setText(f"Exported {len(records)} transcript sequence(s) to {path}")

    def _export_proteins(self) -> None:
        if not self._run_result:
            return
        records = derive_protein_sequences(self._run_result, self._record.sequence)
        if not records:
            self._status.setText("No CDS/protein features available for export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Protein FASTA",
            f"{self._record.id}.proteins.fasta",
            "FASTA (*.fasta *.fa *.faa);;All Files (*)",
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(export_fasta(records))
        self._status.setText(f"Exported {len(records)} protein sequence(s) to {path}")

    def _export_gff3(self) -> None:
        if not self._run_result:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Gene Predictions as GFF3",
            f"{self._record.id}.gff3",
            "GFF3 (*.gff3);;All Files (*)",
        )
        if not path:
            return
        text = export_predictions_gff3(self._run_result, self._record.id)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        self._status.setText(f"Exported GFF3 to {path}")

    def _import_gff3(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Gene Predictions from GFF3",
            "",
            "GFF3 (*.gff3 *.gff);;All Files (*)",
        )
        if not path:
            return

        try:
            with open(path, "r", encoding="utf-8") as handle:
                text = handle.read()
            seq_id, run = import_predictions_gff3(text, preferred_seq_id=self._record.id)
        except Exception as exc:
            self._status.setText(f"GFF3 import failed: {exc}")
            return

        self._run_result = run
        self._selected_index = 0 if run.genes else None
        self._status.setText(
            f"{len(run.genes)} prediction(s) imported from {path} for seqid '{seq_id}'"
        )
        self._refresh_visible_features()

    def _on_row_selected(self) -> None:
        if not self._run_result:
            return

        rows = {self._table.row(i) for i in self._table.selectedItems()}
        if not rows:
            return

        row = min(rows)
        if row >= len(self._visible_indices):
            return

        self._selected_index = self._visible_indices[row]
        gene: GenePrediction = self._run_result.genes[self._selected_index]
        aa_len = len(gene.protein)
        lines = [
            f"Feature type: {gene.feature_type}",
            f"Coordinates: {gene.start + 1}-{gene.end} ({gene.length_nt} nt)",
            f"Strand: {gene.strand}",
            f"Source: {gene.source}",
            f"Feature ID: {gene.feature_id}",
            f"Parent ID: {gene.parent_id}",
            f"Score: {'' if gene.confidence is None else f'{gene.confidence:.3f}'}",
            f"Protein length: {aa_len} aa",
            "",
        ]

        if gene.protein:
            lines.append("Protein translation:")
            lines.extend(gene.protein[i:i + 70] for i in range(0, aa_len, 70))
        else:
            lines.append("Protein translation unavailable for this prediction.")

        self._details.setPlainText("\n".join(lines))
        self._draw_track()

    def _draw_track(self) -> None:
        self._fig.clear()
        ax = self._fig.add_subplot(111)

        seq_len = max(1, len(self._record.sequence.replace("-", "").replace(".", "")))
        ax.set_xlim(0, seq_len)
        ax.set_ylim(-3.3, 3.8)
        ax.set_yticks([3.2, 2.4, 1.6, 0.8, -0.1, -0.9, -1.7, -2.5])
        ax.set_yticklabels([
            "+ gene", "+ mRNA", "+ exon", "+ CDS",
            "- gene", "- mRNA", "- exon", "- CDS",
        ])
        ax.set_xlabel("Position (bp)")
        ax.set_title(f"Gene Feature Track — {self._record.id}", fontweight="bold")
        for y in [3.2, 2.4, 1.6, 0.8, -0.1, -0.9, -1.7, -2.5]:
            ax.axhline(y, color="#E5E5EA", linewidth=0.6, zorder=0)

        genes = [] if not self._run_result else self._run_result.genes
        visible_genes = [genes[idx] for idx in self._visible_indices] if self._run_result else []

        centers: dict[str, tuple[float, float]] = {}
        for idx in self._visible_indices:
            gene = genes[idx]
            y = _FEATURE_Y.get((gene.feature_type, gene.strand), 0.6 if gene.strand == "+" else -2.7)
            h = _FEATURE_HEIGHT.get(gene.feature_type, 0.34)
            face = _TRACK_COLOR.get(gene.strand, "#888888")
            edge = "#1C1C1E" if idx == self._selected_index else "white"
            line_w = 1.4 if idx == self._selected_index else 0.6

            rect = mpatches.FancyBboxPatch(
                (gene.start, y),
                max(1, gene.end - gene.start),
                h,
                boxstyle="round,pad=0,rounding_size=2",
                facecolor=face,
                edgecolor=edge,
                linewidth=line_w,
                alpha=0.9,
            )
            ax.add_patch(rect)
            if gene.feature_id:
                centers[gene.feature_id] = (gene.start + max(1, gene.end - gene.start) / 2, y + h / 2)

            if gene.feature_type in {"gene", "mRNA", "CDS"}:
                arrow_x = gene.end - 8 if gene.strand == "+" else gene.start + 8
                dx = 6 if gene.strand == "+" else -6
                ax.arrow(
                    arrow_x,
                    y + h / 2,
                    dx,
                    0,
                    width=0.01,
                    head_width=0.14,
                    head_length=6,
                    length_includes_head=True,
                    color="#1C1C1E",
                    zorder=3,
                )

        for gene in visible_genes:
            if not gene.parent_id or gene.parent_id not in centers or not gene.feature_id:
                continue
            child_center = centers.get(gene.feature_id)
            parent_center = centers.get(gene.parent_id)
            if not child_center or not parent_center:
                continue

            px, py = parent_center
            cx, cy = child_center
            mid_x = cx
            ax.plot([px, mid_x], [py, py], color="#636366", linewidth=0.9, alpha=0.9, zorder=1)
            ax.plot([mid_x, mid_x], [py, cy], color="#636366", linewidth=0.9, alpha=0.9, zorder=1)

        ax.spines[["top", "right", "left"]].set_visible(False)
        self._canvas.draw()
