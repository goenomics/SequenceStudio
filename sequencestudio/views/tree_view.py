"""Phylogenetic tree view — NJ/UPGMA with full display customisation."""

from __future__ import annotations

from math import log
from typing import Optional

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QLabel, QComboBox, QPushButton, QCheckBox,
    QDoubleSpinBox, QSpinBox, QGroupBox, QFormLayout,
    QFileDialog, QMessageBox, QScrollArea, QSizePolicy,
)
from PySide6.QtCore import Qt

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from sequencestudio.models.sequence import SequenceAlignment


# ---------------------------------------------------------------------------
# Distance functions
# ---------------------------------------------------------------------------

def _p_distance(s1: str, s2: str) -> float:
    """Raw proportion of differing sites (ignoring gaps)."""
    diffs = valid = 0
    for a, b in zip(s1, s2):
        if a == "-" or b == "-":
            continue
        valid += 1
        if a.upper() != b.upper():
            diffs += 1
    return (diffs / valid) if valid else 1.0


def _jukes_cantor(s1: str, s2: str) -> float:
    """Jukes-Cantor corrected distance."""
    p = _p_distance(s1, s2)
    p = min(p, 0.749)   # avoid log(0)
    try:
        return -0.75 * log(1.0 - (4.0 / 3.0) * p)
    except ValueError:
        return 1.0


def _kimura2p(s1: str, s2: str) -> float:
    """Kimura 2-parameter corrected distance (DNA only)."""
    transitions    = {"AG", "GA", "CT", "TC"}
    transversions  = {"AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"}
    ti = tv = valid = 0
    for a, b in zip(s1, s2):
        if a == "-" or b == "-":
            continue
        valid += 1
        pair = a.upper() + b.upper()
        if pair in transitions:
            ti += 1
        elif pair in transversions:
            tv += 1
    if valid == 0:
        return 1.0
    P = ti / valid   # transition proportion
    Q = tv / valid   # transversion proportion
    try:
        return -0.5 * log(1.0 - 2 * P - Q) - 0.25 * log(1.0 - 2 * Q)
    except ValueError:
        return 1.0


_DIST_FUNCS = {
    "p-distance":     _p_distance,
    "Jukes-Cantor":   _jukes_cantor,
    "Kimura 2P":      _kimura2p,
}

_TIP_PALETTES = ["Black", "By clade", "Rainbow", "Tab20"]


def _build_distance_matrix(aln: SequenceAlignment, dist_fn) -> DistanceMatrix:
    names = [r.id for r in aln.records]
    seqs  = [r.sequence for r in aln.records]
    n = len(seqs)
    mat: list[list[float]] = []
    for i in range(n):
        row: list[float] = []
        for j in range(i + 1):
            row.append(0.0 if i == j else dist_fn(seqs[i], seqs[j]))
        mat.append(row)
    return DistanceMatrix(names, mat)


def _tip_label_colors(tree, palette: str) -> dict[str, str]:
    """Build a {name: color_hex} dict for terminal nodes."""
    terminals = [c.name for c in tree.get_terminals() if c.name]
    n = len(terminals)
    if palette == "Black" or n == 0:
        return {}
    cmap_name = {"By clade": "Set2", "Rainbow": "hsv", "Tab20": "tab20"}.get(palette, "tab20")
    cmap = matplotlib.colormaps.get_cmap(cmap_name)
    return {name: matplotlib.colors.to_hex(cmap(i / max(n - 1, 1)))
            for i, name in enumerate(terminals)}


# ---------------------------------------------------------------------------
# Main widget
# ---------------------------------------------------------------------------

class TreeView(QWidget):
    """Phylogenetic tree with a comprehensive options panel."""

    def __init__(self, alignment: SequenceAlignment, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._aln  = alignment
        self._tree = None
        self._build_ui()
        self._draw()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: options panel ────────────────────────────────────────
        opts_outer = QWidget()
        opts_outer.setFixedWidth(230)
        opts_layout = QVBoxLayout(opts_outer)
        opts_layout.setContentsMargins(6, 6, 6, 6)
        opts_layout.setSpacing(8)

        # Construction
        grp_build = QGroupBox("Construction")
        fl = QFormLayout(grp_build)
        fl.setSpacing(5)

        self._method_cb = QComboBox()
        self._method_cb.addItems(["Neighbor-Joining (NJ)", "UPGMA"])
        fl.addRow("Method:", self._method_cb)

        self._dist_cb = QComboBox()
        self._dist_cb.addItems(list(_DIST_FUNCS.keys()))
        self._dist_cb.setToolTip(
            "p-distance: raw substitution fraction\n"
            "Jukes-Cantor: corrects for multiple hits\n"
            "Kimura 2P: separates transitions/transversions (DNA)"
        )
        fl.addRow("Distance:", self._dist_cb)

        opts_layout.addWidget(grp_build)

        # Layout
        grp_layout = QGroupBox("Layout")
        fl2 = QFormLayout(grp_layout)
        fl2.setSpacing(5)

        self._orient_cb = QComboBox()
        self._orient_cb.addItems(["Left → Right", "Right → Left", "Top → Bottom"])
        fl2.addRow("Orientation:", self._orient_cb)

        self._lw_spin = QDoubleSpinBox()
        self._lw_spin.setRange(0.5, 5.0)
        self._lw_spin.setSingleStep(0.25)
        self._lw_spin.setValue(1.5)
        fl2.addRow("Line width:", self._lw_spin)

        self._fontsize_spin = QSpinBox()
        self._fontsize_spin.setRange(5, 18)
        self._fontsize_spin.setValue(9)
        fl2.addRow("Label font:", self._fontsize_spin)

        self._branch_len_cb = QCheckBox("Show branch lengths")
        self._branch_len_cb.setChecked(True)
        fl2.addRow(self._branch_len_cb)

        self._confidence_cb = QCheckBox("Show confidence")
        self._confidence_cb.setChecked(False)
        fl2.addRow(self._confidence_cb)

        opts_layout.addWidget(grp_layout)

        # Colours
        grp_col = QGroupBox("Colours")
        fl3 = QFormLayout(grp_col)
        fl3.setSpacing(5)

        self._tip_col_cb = QComboBox()
        self._tip_col_cb.addItems(_TIP_PALETTES)
        fl3.addRow("Tip labels:", self._tip_col_cb)

        self._branch_col_cb = QComboBox()
        self._branch_col_cb.addItems(["Black", "Steel blue", "Forest green", "Crimson"])
        fl3.addRow("Branches:", self._branch_col_cb)

        opts_layout.addWidget(grp_col)

        # Draw button
        draw_btn = QPushButton("Draw Tree")
        draw_btn.setStyleSheet(
            "QPushButton{background:#007AFF;color:white;border-radius:6px;padding:5px;}"
            "QPushButton:hover{background:#005ecb;}"
        )
        draw_btn.clicked.connect(self._draw)
        opts_layout.addWidget(draw_btn)

        # Export
        grp_export = QGroupBox("Export")
        exp_layout = QVBoxLayout(grp_export)
        exp_layout.setSpacing(4)
        for label, fmt in [("Save PNG…", "png"), ("Save SVG…", "svg"), ("Save PDF…", "pdf")]:
            btn = QPushButton(label)
            btn.clicked.connect(lambda checked, f=fmt: self._export(f))
            exp_layout.addWidget(btn)
        opts_layout.addWidget(grp_export)

        opts_layout.addStretch()
        splitter.addWidget(opts_outer)

        # ── Right: matplotlib canvas ───────────────────────────────────
        canvas_w = QWidget()
        canvas_layout = QVBoxLayout(canvas_w)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        self._fig    = Figure(tight_layout=True)
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        canvas_layout.addWidget(self._canvas)
        splitter.addWidget(canvas_w)

        splitter.setSizes([230, 770])
        root.addWidget(splitter)

    # ------------------------------------------------------------------
    # Tree drawing
    # ------------------------------------------------------------------

    def _draw(self) -> None:
        if self._aln.num_sequences < 3:
            QMessageBox.warning(self, "Tree", "Need at least 3 sequences.")
            return

        # Build tree
        dist_fn = _DIST_FUNCS[self._dist_cb.currentText()]
        try:
            dm          = _build_distance_matrix(self._aln, dist_fn)
            constructor = DistanceTreeConstructor()
            if "UPGMA" in self._method_cb.currentText():
                self._tree = constructor.upgma(dm)
            else:
                self._tree = constructor.nj(dm)
        except Exception as exc:
            QMessageBox.critical(self, "Tree Error", str(exc))
            return

        # Branch labels
        if self._branch_len_cb.isChecked():
            branch_labels = lambda c: (
                f"{c.branch_length:.4f}" if c.branch_length and c.branch_length > 0 else ""
            )
        else:
            branch_labels = lambda c: ""

        # Tip label colours
        tip_colors = _tip_label_colors(self._tree, self._tip_col_cb.currentText())

        # Branch colour
        branch_color_map = {
            "Black": "#000000",
            "Steel blue": "#4682B4",
            "Forest green": "#228B22",
            "Crimson": "#DC143C",
        }
        branch_color = branch_color_map.get(self._branch_col_cb.currentText(), "#000000")

        # Draw
        self._fig.clear()
        ax = self._fig.add_subplot(111)

        Phylo.draw(
            self._tree,
            axes=ax,
            do_show=False,
            label_func=lambda c: c.name if c.name else "",
            label_colors=tip_colors,
            branch_labels=branch_labels,
            show_confidence=self._confidence_cb.isChecked(),
        )

        # Apply matplotlib-level styling
        lw = self._lw_spin.value()
        fs = self._fontsize_spin.value()
        for line in ax.get_lines():
            line.set_linewidth(lw)
            line.set_color(branch_color)
        for txt in ax.texts:
            txt.set_fontsize(fs)
        for lbl in ax.get_yticklabels():
            lbl.set_fontsize(fs)

        # Orientation
        orient = self._orient_cb.currentText()
        if orient == "Right → Left":
            ax.invert_xaxis()
        elif orient == "Top → Bottom":
            ax.invert_yaxis()

        ax.set_title(
            f"{self._aln.name}  —  {self._method_cb.currentText()}  "
            f"({self._dist_cb.currentText()})",
            fontsize=10, fontweight="bold",
        )
        ax.set_xlabel("Branch length")
        ax.spines[["top", "right"]].set_visible(False)

        self._canvas.draw()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export(self, fmt: str) -> None:
        if self._tree is None:
            QMessageBox.information(self, "Export", "Draw the tree first.")
            return
        ext_filter = {"png": "PNG Image (*.png)", "svg": "SVG Vector (*.svg)",
                      "pdf": "PDF Document (*.pdf)"}
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Tree", f"tree.{fmt}", ext_filter.get(fmt, "*.*")
        )
        if not path:
            return
        try:
            self._fig.savefig(path, dpi=200, bbox_inches="tight")
            QMessageBox.information(self, "Export", f"Saved to {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export Error", str(exc))

