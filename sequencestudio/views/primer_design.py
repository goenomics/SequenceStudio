"""Primer design dialog."""

from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QLabel, QSpinBox, QDoubleSpinBox, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QTextEdit,
    QSplitter, QMessageBox,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from sequencestudio.models.sequence import SequenceRecord
from sequencestudio.core.primer import design_primers, calculate_tm, gc_content, Primer


class PrimerDesignDialog(QDialog):
    """Interactive primer design with nearest-neighbor Tm calculation."""

    def __init__(self, record: SequenceRecord, parent=None) -> None:
        super().__init__(parent)
        self._record = record
        self._primers: list[Primer] = []
        self.setWindowTitle(f"Primer Design — {record.id}")
        self.setMinimumSize(900, 620)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # ---- Parameters ----
        params = QGroupBox("Parameters")
        form = QFormLayout(params)

        seq_len = len(self._record.sequence.replace("-", ""))

        self._tgt_start = QSpinBox()
        self._tgt_start.setRange(1, seq_len)
        self._tgt_start.setValue(max(1, seq_len // 3))
        self._tgt_start.setSuffix(" bp")
        form.addRow("Target start (1-based):", self._tgt_start)

        self._tgt_end = QSpinBox()
        self._tgt_end.setRange(1, seq_len)
        self._tgt_end.setValue(min(seq_len, seq_len * 2 // 3))
        self._tgt_end.setSuffix(" bp")
        form.addRow("Target end (1-based):", self._tgt_end)

        len_row = QHBoxLayout()
        self._min_len = QSpinBox()
        self._min_len.setRange(15, 35)
        self._min_len.setValue(18)
        self._max_len = QSpinBox()
        self._max_len.setRange(15, 40)
        self._max_len.setValue(25)
        len_row.addWidget(self._min_len)
        len_row.addWidget(QLabel(" – "))
        len_row.addWidget(self._max_len)
        len_row.addStretch()
        form.addRow("Primer length (nt):", len_row)

        tm_row = QHBoxLayout()
        self._tm_min = QDoubleSpinBox()
        self._tm_min.setRange(30.0, 90.0)
        self._tm_min.setValue(55.0)
        self._tm_min.setSuffix(" °C")
        self._tm_max = QDoubleSpinBox()
        self._tm_max.setRange(30.0, 90.0)
        self._tm_max.setValue(65.0)
        self._tm_max.setSuffix(" °C")
        tm_row.addWidget(self._tm_min)
        tm_row.addWidget(QLabel(" – "))
        tm_row.addWidget(self._tm_max)
        tm_row.addStretch()
        form.addRow("Tm range:", tm_row)

        gc_row = QHBoxLayout()
        self._gc_min = QDoubleSpinBox()
        self._gc_min.setRange(0, 100)
        self._gc_min.setValue(40.0)
        self._gc_min.setSuffix("%")
        self._gc_max = QDoubleSpinBox()
        self._gc_max.setRange(0, 100)
        self._gc_max.setValue(60.0)
        self._gc_max.setSuffix("%")
        gc_row.addWidget(self._gc_min)
        gc_row.addWidget(QLabel(" – "))
        gc_row.addWidget(self._gc_max)
        gc_row.addStretch()
        form.addRow("GC% range:", gc_row)

        root.addWidget(params)

        btn_row = QHBoxLayout()
        run_btn = QPushButton("Design Primers")
        run_btn.setFixedHeight(32)
        run_btn.setStyleSheet(
            "QPushButton{background:#34C759;color:white;border-radius:6px;padding:0 20px;}"
            "QPushButton:hover{background:#248a3d;}"
        )
        run_btn.clicked.connect(self._run)
        btn_row.addStretch()
        btn_row.addWidget(run_btn)
        root.addLayout(btn_row)

        self._status = QLabel("Set target region and click Design Primers.")
        self._status.setObjectName("status_label")
        root.addWidget(self._status)

        # ---- Results ----
        sp = QSplitter(Qt.Orientation.Vertical)

        self._table = QTableWidget(0, 6)
        self._table.setHorizontalHeaderLabels(
            ["Name", "Strand", "Sequence (5'→3')", "Length", "Tm (°C)", "GC%"]
        )
        self._table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.itemSelectionChanged.connect(self._show_detail)
        sp.addWidget(self._table)

        self._detail = QTextEdit()
        self._detail.setReadOnly(True)
        self._detail.setFont(QFont("Menlo", 10))
        self._detail.setPlaceholderText("Select a primer to view details.")
        sp.addWidget(self._detail)
        sp.setSizes([320, 160])
        root.addWidget(sp)

    def _run(self) -> None:
        ts = self._tgt_start.value() - 1   # convert to 0-based
        te = self._tgt_end.value()
        if ts >= te:
            QMessageBox.warning(self, "Invalid Range", "Target start must be before target end.")
            return

        template = self._record.sequence.replace("-", "")
        try:
            self._primers = design_primers(
                template, ts, te,
                length_range=(self._min_len.value(), self._max_len.value()),
                tm_range=(self._tm_min.value(), self._tm_max.value()),
                gc_range=(self._gc_min.value(), self._gc_max.value()),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Error", str(exc))
            return

        self._table.setRowCount(len(self._primers))
        for row, p in enumerate(self._primers):
            strand_lbl = "Forward (+)" if p.strand == 1 else "Reverse (−)"
            self._table.setItem(row, 0, QTableWidgetItem(p.name))
            self._table.setItem(row, 1, QTableWidgetItem(strand_lbl))
            self._table.setItem(row, 2, QTableWidgetItem(p.sequence))
            self._table.setItem(row, 3, QTableWidgetItem(str(p.length)))
            self._table.setItem(row, 4, QTableWidgetItem(str(p.tm)))
            self._table.setItem(row, 5, QTableWidgetItem(f"{p.gc_pct:.1f}"))

        n_fwd = sum(1 for p in self._primers if p.strand == 1)
        n_rev = len(self._primers) - n_fwd
        self._status.setText(f"Found {n_fwd} forward + {n_rev} reverse primers.")

    def _show_detail(self) -> None:
        rows = {self._table.row(i) for i in self._table.selectedItems()}
        if not rows:
            return
        p = self._primers[min(rows)]
        strand_lbl = "Forward (+)" if p.strand == 1 else "Reverse (−)"
        self._detail.setPlainText(
            f"Name:     {p.name}\n"
            f"Strand:   {strand_lbl}\n"
            f"Sequence: 5'-{p.sequence}-3'\n"
            f"Length:   {p.length} nt\n"
            f"Tm:       {p.tm} °C\n"
            f"GC:       {p.gc_pct:.1f}%\n"
            f"Position: {p.start + 1}–{p.end} (on template forward strand)\n"
        )
