"""Pairwise dot-plot view using matplotlib."""

from __future__ import annotations

from typing import Optional

import numpy as np
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QPushButton, QProgressBar,
)
from PySide6.QtCore import Qt, QThread, Signal

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from sequencestudio.models.sequence import SequenceRecord


def _compute_dot(seq1: str, seq2: str, window: int) -> np.ndarray:
    """Compute a dot-plot matrix using a sliding window."""
    s1 = seq1.upper().replace("-", "").replace(".", "")
    s2 = seq2.upper().replace("-", "").replace(".", "")
    n, m = len(s1), len(s2)
    mat = np.zeros((n, m), dtype=np.float32)
    half = window // 2
    threshold = window // 2 + 1  # minimum matches to mark a dot

    for i in range(half, n - half):
        for j in range(half, m - half):
            matches = sum(
                s1[i - half + k] == s2[j - half + k]
                for k in range(window)
            )
            if matches >= threshold:
                mat[i, j] = 1.0
    return mat


class _DotWorker(QThread):
    finished = Signal(object, str, str)   # ndarray, name1, name2
    error    = Signal(str)

    def __init__(self, s1: str, n1: str, s2: str, n2: str, window: int) -> None:
        super().__init__()
        self._s1, self._n1 = s1, n1
        self._s2, self._n2 = s2, n2
        self._window = window

    def run(self) -> None:
        try:
            mat = _compute_dot(self._s1, self._s2, self._window)
            self.finished.emit(mat, self._n1, self._n2)
        except Exception as exc:
            self.error.emit(str(exc))


class DotPlotView(QWidget):
    """Dot-plot comparing two sequences from a loaded set."""

    def __init__(self, records: list[SequenceRecord], parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._records = records
        self._worker: Optional[_DotWorker] = None
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # Sequence selectors
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("X-axis:"))
        self._seq1_cb = QComboBox()
        for r in self._records:
            self._seq1_cb.addItem(r.id)
        ctrl.addWidget(self._seq1_cb)

        ctrl.addWidget(QLabel("  Y-axis:"))
        self._seq2_cb = QComboBox()
        for r in self._records:
            self._seq2_cb.addItem(r.id)
        if len(self._records) > 1:
            self._seq2_cb.setCurrentIndex(1)
        ctrl.addWidget(self._seq2_cb)

        ctrl.addWidget(QLabel("  Window:"))
        self._win_spin = QSpinBox()
        self._win_spin.setRange(1, 51)
        self._win_spin.setSingleStep(2)
        self._win_spin.setValue(7)
        ctrl.addWidget(self._win_spin)

        self._run_btn = QPushButton("Draw Dot Plot")
        self._run_btn.setStyleSheet(
            "QPushButton{background:#007AFF;color:white;border-radius:6px;padding:3px 14px;}"
            "QPushButton:hover{background:#005ecb;}"
            "QPushButton:disabled{background:#AEAEB2;}"
        )
        self._run_btn.clicked.connect(self._run)
        ctrl.addWidget(self._run_btn)
        ctrl.addStretch()
        root.addLayout(ctrl)

        self._progress = QProgressBar()
        self._progress.setVisible(False)
        root.addWidget(self._progress)

        self._status = QLabel(
            "Select two sequences and click 'Draw Dot Plot'.\n"
            "Large sequences may take a moment."
        )
        self._status.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._status.setObjectName("status_label")
        root.addWidget(self._status)

        # Matplotlib canvas
        self._fig = Figure(tight_layout=True)
        self._canvas = FigureCanvas(self._fig)
        root.addWidget(self._canvas)

    def _run(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        i1 = self._seq1_cb.currentIndex()
        i2 = self._seq2_cb.currentIndex()
        r1, r2 = self._records[i1], self._records[i2]
        window = self._win_spin.value()

        # Warn about very large sequences
        max_len = 5000
        if max(r1.length, r2.length) > max_len:
            self._status.setText(
                f"Warning: sequences > {max_len:,} bp may be slow. Drawing anyway…"
            )

        self._run_btn.setEnabled(False)
        self._progress.setVisible(True)
        self._progress.setRange(0, 0)

        self._worker = _DotWorker(r1.sequence, r1.id, r2.sequence, r2.id, window)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, mat: np.ndarray, name1: str, name2: str) -> None:
        self._run_btn.setEnabled(True)
        self._progress.setVisible(False)

        self._fig.clear()
        ax = self._fig.add_subplot(111)
        ax.imshow(mat.T, origin="lower", aspect="auto",
                  cmap="binary", interpolation="nearest")
        ax.set_xlabel(f"{name1} →")
        ax.set_ylabel(f"{name2} →")
        ax.set_title(f"Dot Plot  |  window = {self._win_spin.value()}")
        self._canvas.draw()
        self._status.setText(
            f"Dot plot: {name1} ({mat.shape[0]} bp) vs {name2} ({mat.shape[1]} bp)"
        )

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._progress.setVisible(False)
        self._status.setText(f"Error: {msg}")
