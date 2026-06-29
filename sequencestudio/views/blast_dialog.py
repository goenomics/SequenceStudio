"""NCBI BLAST search dialog."""

from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QLabel, QComboBox, QSpinBox, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QTextEdit, QSplitter, QMessageBox,
)
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QFont

from sequencestudio.models.sequence import SequenceRecord
from sequencestudio.core.blast import run_blast_ncbi, BlastResult, BLAST_PROGRAMS, BLAST_DATABASES


class _BlastWorker(QThread):
    finished = Signal(object)   # BlastResult
    progress = Signal(str)
    error    = Signal(str)

    def __init__(self, seq: str, sid: str, prog: str, db: str, n: int) -> None:
        super().__init__()
        self._seq, self._sid = seq, sid
        self._prog, self._db, self._n = prog, db, n

    def run(self) -> None:
        try:
            res = run_blast_ncbi(
                self._seq, self._sid, self._prog, self._db, self._n,
                self.progress.emit,
            )
            self.finished.emit(res)
        except Exception as exc:
            self.error.emit(str(exc))


class BlastDialog(QDialog):
    """NCBI remote BLAST search interface."""

    def __init__(self, record: SequenceRecord, parent=None) -> None:
        super().__init__(parent)
        self._record = record
        self._result: Optional[BlastResult] = None
        self._worker: Optional[_BlastWorker] = None
        self.setWindowTitle(f"BLAST — {record.id}")
        self.setMinimumSize(940, 680)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # Parameters
        params = QGroupBox("Search Parameters")
        form  = QFormLayout(params)

        self._prog_cb = QComboBox()
        self._prog_cb.addItems(BLAST_PROGRAMS)
        default_prog = "blastp" if self._record.molecule_type == "Protein" else "blastn"
        self._prog_cb.setCurrentText(default_prog)
        self._prog_cb.currentTextChanged.connect(self._update_dbs)
        form.addRow("Program:", self._prog_cb)

        self._db_cb = QComboBox()
        form.addRow("Database:", self._db_cb)
        self._update_dbs(default_prog)

        self._hits_spin = QSpinBox()
        self._hits_spin.setRange(1, 500)
        self._hits_spin.setValue(50)
        form.addRow("Max hits:", self._hits_spin)

        mol = self._record.molecule_type
        form.addRow("Query:",
                    QLabel(f"{self._record.id}  ({self._record.length:,} bp, {mol})"))
        root.addWidget(params)

        # Run button
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton("Run BLAST (NCBI)")
        self._run_btn.setFixedHeight(32)
        self._run_btn.setStyleSheet(
            "QPushButton{background:#007AFF;color:white;border-radius:6px;padding:0 20px;}"
            "QPushButton:hover{background:#005ecb;}"
            "QPushButton:disabled{background:#AEAEB2;}"
        )
        self._run_btn.clicked.connect(self._run)
        btn_row.addStretch()
        btn_row.addWidget(self._run_btn)
        root.addLayout(btn_row)

        self._status = QLabel("Ready — click 'Run BLAST' to query NCBI (may take 30–90 s).")
        self._status.setObjectName("status_label")
        root.addWidget(self._status)

        # Results
        sp = QSplitter(Qt.Orientation.Vertical)

        self._table = QTableWidget(0, 7)
        self._table.setHorizontalHeaderLabels(
            ["Description", "Accession", "Length", "Score", "E-value", "Identity", "Ident %"]
        )
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.itemSelectionChanged.connect(self._show_detail)
        sp.addWidget(self._table)

        self._detail = QTextEdit()
        self._detail.setReadOnly(True)
        self._detail.setFont(QFont("Menlo", 10))
        self._detail.setPlaceholderText("Select a hit to view alignment.")
        sp.addWidget(self._detail)

        sp.setSizes([320, 200])
        root.addWidget(sp)

    def _update_dbs(self, prog: str) -> None:
        self._db_cb.clear()
        self._db_cb.addItems(BLAST_DATABASES.get(prog, ["nt"]))

    def _run(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        self._run_btn.setEnabled(False)
        self._table.setRowCount(0)
        self._detail.clear()
        self._status.setText("Submitting to NCBI…")

        self._worker = _BlastWorker(
            self._record.sequence, self._record.id,
            self._prog_cb.currentText(), self._db_cb.currentText(),
            self._hits_spin.value(),
        )
        self._worker.finished.connect(self._on_done)
        self._worker.progress.connect(self._status.setText)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, result: BlastResult) -> None:
        self._run_btn.setEnabled(True)
        self._result = result
        self._table.setRowCount(len(result.hits))
        for row, hit in enumerate(result.hits):
            self._table.setItem(row, 0, QTableWidgetItem(hit.title))
            self._table.setItem(row, 1, QTableWidgetItem(hit.accession))
            self._table.setItem(row, 2, QTableWidgetItem(f"{hit.length:,}"))
            self._table.setItem(row, 3, QTableWidgetItem(str(int(hit.score))))
            self._table.setItem(row, 4, QTableWidgetItem(f"{hit.evalue:.2e}"))
            self._table.setItem(row, 5, QTableWidgetItem(str(hit.identity)))
            self._table.setItem(row, 6, QTableWidgetItem(f"{hit.identity_pct:.1f}%"))
        self._status.setText(f"{result.num_hits} hit(s) found.")

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status.setText(f"Error: {msg}")
        QMessageBox.critical(self, "BLAST Error", msg)

    def _show_detail(self) -> None:
        if not self._result:
            return
        rows = {self._table.row(i) for i in self._table.selectedItems()}
        if not rows:
            return
        hit = self._result.hits[min(rows)]
        self._detail.setPlainText(
            f">{hit.title}\n"
            f"Accession: {hit.accession}   Length: {hit.length:,}\n"
            f"Score: {hit.score}   E-value: {hit.evalue:.2e}   "
            f"Identity: {hit.identity} ({hit.identity_pct:.1f}%)\n\n"
            f"{hit.alignment_text}"
        )
