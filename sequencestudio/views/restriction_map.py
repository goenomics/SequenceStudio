"""Restriction enzyme map view using matplotlib."""

from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QCheckBox,
    QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
    QHeaderView,
)
from PySide6.QtCore import Qt

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from sequencestudio.models.sequence import SequenceRecord
from sequencestudio.core.restriction import analyze_sequence, RestrictionResult


class RestrictionMapView(QWidget):
    """Visual restriction-enzyme map for a single sequence."""

    def __init__(self, record: SequenceRecord, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._record = record
        self._results: list[RestrictionResult] = []
        self._selected: set[str] = set()
        self._build_ui()
        self._analyze()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # Controls bar
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel(
            f"<b>{self._record.id}</b>  ({self._record.length:,} bp)"
        ))
        ctrl.addStretch()
        self._linear_cb = QCheckBox("Linear topology")
        self._linear_cb.setChecked(True)
        self._linear_cb.toggled.connect(self._analyze)
        ctrl.addWidget(self._linear_cb)
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self._analyze)
        ctrl.addWidget(refresh)
        root.addLayout(ctrl)

        # Splitter: table | map
        sp = QSplitter(Qt.Orientation.Horizontal)

        # Results table
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(QLabel("<b>Restriction Sites</b>"))
        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(["Enzyme", "Cuts", "Fragments (bp)"])
        self._table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QTableWidget.SelectionMode.MultiSelection)
        self._table.itemSelectionChanged.connect(self._on_enzyme_select)
        ll.addWidget(self._table)
        left.setMinimumWidth(260)
        sp.addWidget(left)

        # Matplotlib canvas
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        self._fig = Figure(figsize=(9, 3.5), tight_layout=True)
        self._canvas = FigureCanvas(self._fig)
        rl.addWidget(self._canvas)
        sp.addWidget(right)

        sp.setSizes([270, 730])
        root.addWidget(sp)

    def _analyze(self) -> None:
        linear = self._linear_cb.isChecked()
        try:
            self._results = analyze_sequence(
                self._record.sequence, linear=linear
            )
        except Exception:
            self._results = []
        self._populate_table()
        self._selected.clear()
        self._draw_map()

    def _populate_table(self) -> None:
        self._table.setRowCount(len(self._results))
        for row, r in enumerate(self._results):
            self._table.setItem(row, 0, QTableWidgetItem(r.enzyme))
            self._table.setItem(row, 1, QTableWidgetItem(str(r.num_cuts)))
            frags = ", ".join(str(f) for f in r.fragment_sizes[:6])
            if len(r.fragment_sizes) > 6:
                frags += " …"
            self._table.setItem(row, 2, QTableWidgetItem(frags))

    def _on_enzyme_select(self) -> None:
        rows = {self._table.row(i) for i in self._table.selectedItems()}
        self._selected = {
            self._results[r].enzyme for r in rows if r < len(self._results)
        }
        self._draw_map()

    def _draw_map(self) -> None:
        self._fig.clear()
        ax = self._fig.add_subplot(111)
        seq_len = len(self._record.sequence.replace("-", ""))

        # Backbone
        ax.axhline(0.5, color="#1C1C1E", linewidth=4, zorder=1)
        ax.set_xlim(0, seq_len)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("Position (bp)")
        ax.set_title(f"Restriction Map — {self._record.id}", fontweight="bold")
        ax.yaxis.set_visible(False)
        ax.spines[["top", "left", "right"]].set_visible(False)

        if self._selected:
            to_show = [r for r in self._results if r.enzyme in self._selected]
        else:
            to_show = [r for r in self._results if r.num_cuts <= 3][:12]

        cmap = matplotlib.colormaps.get_cmap("tab20")
        for i, res in enumerate(to_show):
            color = cmap(i % 20)
            for pos in res.cut_positions:
                ax.axvline(pos, color=color, linewidth=1.5, alpha=0.85, zorder=2)
                voffset = 0.82 - (i % 6) * 0.07
                ax.text(pos, voffset, res.enzyme,
                        fontsize=6, rotation=90, ha="center", va="bottom",
                        color=color, fontweight="bold")

        self._canvas.draw()
