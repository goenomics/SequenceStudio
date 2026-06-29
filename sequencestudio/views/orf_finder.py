"""ORF Finder — graphical 6-frame map + table."""

from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox,
    QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
    QHeaderView, QTextEdit,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as mpatches

from sequencestudio.models.sequence import SequenceRecord
from sequencestudio.core.orf import find_orfs, ORF

_FRAME_COLOR = {
    1:  "#4C9BE8",  # forward frames — blue shades
    2:  "#34AADC",
    3:  "#5AC8FA",
    -1: "#FF6B6B",  # reverse frames — red shades
    -2: "#FF3A30",
    -3: "#FF9500",
}
_FRAME_Y = {
    1: 5, 2: 4, 3: 3, -1: 2, -2: 1, -3: 0
}


class OrfFinderView(QWidget):
    """Six-frame ORF map with sortable result table."""

    def __init__(self, record: SequenceRecord, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._record = record
        self._orfs: list[ORF] = []
        self._build_ui()
        self._run()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # Controls
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel(
            f"<b>{self._record.id}</b>  ({self._record.length:,} bp)"
        ))
        ctrl.addStretch()
        ctrl.addWidget(QLabel("Min ORF length (nt):"))
        self._min_spin = QSpinBox()
        self._min_spin.setRange(30, 10000)
        self._min_spin.setSingleStep(30)
        self._min_spin.setValue(150)
        ctrl.addWidget(self._min_spin)
        run_btn = QPushButton("Find ORFs")
        run_btn.clicked.connect(self._run)
        ctrl.addWidget(run_btn)
        root.addLayout(ctrl)

        # Splitter: map top | table+detail bottom
        outer_sp = QSplitter(Qt.Orientation.Vertical)

        # Matplotlib 6-frame map
        map_w = QWidget()
        ml = QVBoxLayout(map_w)
        ml.setContentsMargins(0, 0, 0, 0)
        self._fig = Figure(figsize=(12, 3), tight_layout=True)
        self._canvas = FigureCanvas(self._fig)
        ml.addWidget(self._canvas)
        outer_sp.addWidget(map_w)

        # Table + detail splitter
        inner_sp = QSplitter(Qt.Orientation.Horizontal)

        # ORF table
        self._table = QTableWidget(0, 5)
        self._table.setHorizontalHeaderLabels(
            ["Frame", "Start", "End", "Length (nt)", "AA Length"]
        )
        for col in range(5):
            self._table.horizontalHeader().setSectionResizeMode(
                col, QHeaderView.ResizeMode.ResizeToContents
            )
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.itemSelectionChanged.connect(self._on_row_selected)
        inner_sp.addWidget(self._table)

        # Protein translation view
        self._protein_view = QTextEdit()
        self._protein_view.setReadOnly(True)
        self._protein_view.setFont(QFont("Menlo", 10))
        self._protein_view.setPlaceholderText("Select an ORF to view translation.")
        inner_sp.addWidget(self._protein_view)
        inner_sp.setSizes([400, 500])

        outer_sp.addWidget(inner_sp)
        outer_sp.setSizes([200, 400])
        root.addWidget(outer_sp)

    def _run(self) -> None:
        min_len = self._min_spin.value()
        clean = self._record.sequence.replace("-", "").replace(".", "")
        self._orfs = find_orfs(clean, min_len)
        self._populate_table()
        self._draw_map()

    def _populate_table(self) -> None:
        self._table.setRowCount(len(self._orfs))
        for row, orf in enumerate(self._orfs):
            frame_lbl = f"+{orf.frame}" if orf.frame > 0 else str(orf.frame)
            self._table.setItem(row, 0, QTableWidgetItem(frame_lbl))
            self._table.setItem(row, 1, QTableWidgetItem(str(orf.start + 1)))
            self._table.setItem(row, 2, QTableWidgetItem(str(orf.end)))
            self._table.setItem(row, 3, QTableWidgetItem(str(orf.length)))
            self._table.setItem(row, 4, QTableWidgetItem(str(len(orf.protein))))

    def _on_row_selected(self) -> None:
        rows = {self._table.row(i) for i in self._table.selectedItems()}
        if not rows:
            return
        orf = self._orfs[min(rows)]
        frame_lbl = f"+{orf.frame}" if orf.frame > 0 else str(orf.frame)
        self._protein_view.setPlainText(
            f"Frame {frame_lbl} | Position {orf.start + 1}–{orf.end} | {orf.length} nt | {len(orf.protein)} aa\n\n"
            + "\n".join(orf.protein[i:i+60] for i in range(0, len(orf.protein), 60))
        )

    def _draw_map(self) -> None:
        self._fig.clear()
        ax = self._fig.add_subplot(111)
        seq_len = len(self._record.sequence.replace("-", ""))

        # Frame labels
        frame_labels = {5: "+1", 4: "+2", 3: "+3", 2: "−1", 1: "−2", 0: "−3"}
        for y, lbl in frame_labels.items():
            ax.text(-seq_len * 0.01, y + 0.5, lbl, ha="right", va="center",
                    fontsize=8, fontweight="bold")

        # Guide lines per frame
        for yi in range(6):
            ax.axhline(yi + 0.5, color="#EEEEEE", linewidth=0.5, zorder=0)

        # ORF boxes
        for orf in self._orfs:
            fy = _FRAME_Y.get(orf.frame, 0)
            color = _FRAME_COLOR.get(orf.frame, "#888888")
            rect = mpatches.FancyBboxPatch(
                (orf.start, fy + 0.1),
                orf.end - orf.start, 0.8,
                boxstyle="round,pad=0",
                facecolor=color, edgecolor="white", linewidth=0.5,
                alpha=0.85,
            )
            ax.add_patch(rect)

        ax.set_xlim(-seq_len * 0.02, seq_len * 1.01)
        ax.set_ylim(-0.1, 6.1)
        ax.set_xlabel("Position (bp)")
        ax.set_title(f"6-frame ORF Map — {self._record.id}", fontweight="bold")
        ax.yaxis.set_visible(False)
        ax.spines[["top", "left", "right"]].set_visible(False)

        self._canvas.draw()
