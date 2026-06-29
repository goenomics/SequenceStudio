"""Multiple Sequence Alignment viewer widget."""

from __future__ import annotations

from typing import Optional

from PySide6.QtWidgets import (
    QWidget, QScrollArea, QVBoxLayout, QHBoxLayout,
    QComboBox, QLabel, QSlider, QSizePolicy,
)
from PySide6.QtCore import Qt, QRect, Signal
from PySide6.QtGui import (
    QPainter, QFont, QColor, QPen,
    QMouseEvent, QPaintEvent,
)

from sequencestudio.models.sequence import SequenceAlignment
from sequencestudio.utils.colors import (
    NUCLEOTIDE_COLORS,
    COLOR_SCHEMES, get_text_color,
)


# ---------------------------------------------------------------------------
# Raw drawing canvas
# ---------------------------------------------------------------------------

class _AlignmentCanvas(QWidget):
    """Renders an MSA using only the visible viewport cells each frame."""

    hovered = Signal(int, int)   # row, col (0-based)

    _NAME_W  = 200
    _RULER_H = 22
    _CONS_H  = 22
    _BAR_H   = 28

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._aln: Optional[SequenceAlignment] = None
        self._scheme: dict[str, str] = NUCLEOTIDE_COLORS

        self._cw = 14    # cell width  px
        self._ch = 18    # cell height px

        self._font      = QFont("Menlo", 10, QFont.Weight.Bold)
        self._name_font = QFont("Menlo", max(7, self._ch - 6))

        # Pre-built pens (never recreated)
        self._pen_black = QPen(QColor("#000000"))
        self._pen_white = QPen(QColor("#FFFFFF"))
        self._pen_grey  = QPen(QColor("#DDDDDD"))
        self._pen_dark  = QPen(QColor("#1C1C1E"))
        self._pen_ruler = QPen(QColor("#555555"))

        # Pre-built name-label background colours
        self._bg_even = QColor("#F2F2F7")
        self._bg_odd  = QColor("#FFFFFF")
        self._bg_sel  = QColor("#CCE5FF")

        # Per-scheme colour lookup: char → (bg QColor, True if bg is light)
        # True → use black text pen; False → use white text pen
        self._clut: dict[str, tuple[QColor, bool]] = {}
        self._clut_fallback: tuple[QColor, bool] = (QColor("#DDDDDD"), True)
        self._rebuild_clut()

        # Cached alignment-level data (invalidated on set_alignment / set_scheme)
        self._consensus:    Optional[str]         = None
        self._conservation: Optional[list[float]] = None

        self._sel_col: Optional[int] = None
        self._sel_row: Optional[int] = None
        self.setMouseTracking(True)

    # ------------------------------------------------------------------
    # Public setters
    # ------------------------------------------------------------------

    def set_alignment(self, aln: SequenceAlignment) -> None:
        self._aln = aln
        self._invalidate_cache()
        self._update_size()
        self.update()

    def set_scheme(self, scheme: dict[str, str]) -> None:
        self._scheme = scheme
        self._rebuild_clut()
        self._invalidate_cache()
        self.update()

    def set_cell_size(self, cw: int, ch: int) -> None:
        self._cw, self._ch = cw, ch
        pt = max(6, ch - 5)
        self._font      = QFont("Menlo", pt, QFont.Weight.Bold)
        self._name_font = QFont("Menlo", max(7, ch - 6))
        self._update_size()
        self.update()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _rebuild_clut(self) -> None:
        """Build char → (QColor, is_light_bg) lookup for the current scheme."""
        self._clut = {}
        for char, hex_bg in self._scheme.items():
            bg = QColor(hex_bg)
            light = get_text_color(hex_bg) == "#000000"
            self._clut[char.upper()] = (bg, light)
        self._clut_fallback = (QColor("#DDDDDD"), True)

    def _invalidate_cache(self) -> None:
        self._consensus    = None
        self._conservation = None

    def _ensure_cache(self) -> None:
        if self._aln is None:
            return
        if self._consensus is None:
            self._consensus = self._aln.consensus()
        if self._conservation is None:
            self._conservation = self._aln.conservation()

    def _update_size(self) -> None:
        if not self._aln:
            return
        w = self._NAME_W + self._aln.alignment_length * self._cw + 40
        h = (self._RULER_H
             + self._aln.num_sequences * self._ch
             + self._CONS_H + self._BAR_H + 20)
        self.setMinimumSize(w, h)
        self.resize(w, h)

    def _cell_colors(self, char: str) -> tuple[QColor, bool]:
        return self._clut.get(char.upper(), self._clut_fallback)

    # ------------------------------------------------------------------
    # Paint
    # ------------------------------------------------------------------

    def paintEvent(self, event: QPaintEvent) -> None:
        if not self._aln or not self._aln.records:
            return

        self._ensure_cache()

        p   = QPainter(self)
        clip = event.rect()
        nw, cw, ch = self._NAME_W, self._cw, self._ch
        rh, csh, barh = self._RULER_H, self._CONS_H, self._BAR_H
        alen = self._aln.alignment_length
        nseq = self._aln.num_sequences

        # Don't draw residue text when cells are too small to be legible
        draw_text = (cw >= 8)

        # ── Compute exactly the visible row / column range ──────────────
        col_start = max(0,    (clip.left()   - nw) // cw)
        col_end   = min(alen, (clip.right()  - nw) // cw + 2)
        row_start = max(0,    (clip.top()    - rh) // ch)
        row_end   = min(nseq, (clip.bottom() - rh) // ch + 2)

        # ── Ruler ───────────────────────────────────────────────────────
        if clip.top() < rh:
            p.setFont(QFont("Menlo", 7))
            p.setPen(self._pen_ruler)
            first_tick = (col_start // 10) * 10
            for col in range(first_tick, col_end, 10):
                x = nw + col * cw
                p.drawLine(x, rh - 7, x, rh - 2)
                p.drawText(x + 1, rh - 2, str(col + 1))

        # ── Sequence rows ───────────────────────────────────────────────
        p.setFont(self._font)
        sel_col = self._sel_col

        for ri in range(row_start, row_end):
            rec = self._aln.records[ri]
            yt  = rh + ri * ch

            # Name label (left panel)
            if clip.left() < nw:
                bg_name = (self._bg_sel if self._sel_row == ri
                           else self._bg_even if ri % 2 == 0
                           else self._bg_odd)
                p.fillRect(QRect(0, yt, nw - 2, ch), bg_name)
                p.setPen(self._pen_dark)
                p.setFont(self._name_font)
                p.drawText(QRect(5, yt, nw - 10, ch),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           rec.id[:30])
                p.setFont(self._font)
                p.setPen(self._pen_grey)
                p.drawLine(nw - 1, yt, nw - 1, yt + ch)

            # Residue cells — only the visible column slice
            seq     = rec.sequence
            seq_len = len(seq)

            for ci in range(col_start, col_end):
                char = seq[ci] if ci < seq_len else "-"
                bg, is_light = self._cell_colors(char)

                x = nw + ci * cw

                # Selected-column highlight: compute lighter QColor in-place
                if sel_col == ci:
                    r2, g2, b2 = (
                        min(255, int(bg.red()   * 1.35)),
                        min(255, int(bg.green() * 1.35)),
                        min(255, int(bg.blue()  * 1.35)),
                    )
                    p.fillRect(QRect(x, yt, cw - 1, ch - 1), QColor(r2, g2, b2))
                else:
                    p.fillRect(QRect(x, yt, cw - 1, ch - 1), bg)

                if draw_text:
                    p.setPen(self._pen_black if is_light else self._pen_white)
                    p.drawText(QRect(x, yt, cw, ch),
                               Qt.AlignmentFlag.AlignCenter, char.upper())

        # ── Consensus row ───────────────────────────────────────────────
        cy = rh + nseq * ch
        if clip.bottom() >= cy and self._consensus:
            p.fillRect(QRect(0, cy, nw - 2, csh), QColor("#E8E8ED"))
            p.setFont(QFont("Menlo", max(6, csh - 8), QFont.Weight.Bold))
            p.setPen(self._pen_dark)
            p.drawText(QRect(5, cy, nw - 10, csh),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       "Consensus")
            p.setFont(self._font)
            consensus = self._consensus
            for ci in range(col_start, min(col_end, len(consensus))):
                char = consensus[ci]
                bg, is_light = self._cell_colors(char)
                x = nw + ci * cw
                p.fillRect(QRect(x, cy, cw - 1, csh - 1), bg)
                if draw_text:
                    p.setPen(self._pen_black if is_light else self._pen_white)
                    p.drawText(QRect(x, cy, cw, csh),
                               Qt.AlignmentFlag.AlignCenter, char.upper())

        # ── Conservation bar ────────────────────────────────────────────
        by = cy + csh
        if clip.bottom() >= by and self._conservation:
            p.fillRect(QRect(0, by, nw - 2, barh), QColor("#F2F2F7"))
            p.setFont(QFont("Menlo", 6))
            p.setPen(self._pen_dark)
            p.drawText(QRect(5, by, nw - 10, barh),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       "Conservation")
            conservation = self._conservation
            for ci in range(col_start, min(col_end, len(conservation))):
                score = conservation[ci]
                x     = nw + ci * cw
                bh    = max(1, int(score * (barh - 4)))
                p.fillRect(
                    QRect(x, by + barh - bh - 2, cw - 1, bh),
                    QColor(int(50 + 200 * (1 - score)), int(50 + 200 * score), 80),
                )

        p.end()

    # ------------------------------------------------------------------
    # Mouse
    # ------------------------------------------------------------------

    def mousePressEvent(self, event: QMouseEvent) -> None:
        self._update_selection(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self._update_selection(event, emit_only=True)

    def _update_selection(self, event: QMouseEvent, emit_only: bool = False) -> None:
        if not self._aln:
            return
        x = int(event.position().x())
        y = int(event.position().y())
        col = (x - self._NAME_W) // self._cw
        row = (y - self._RULER_H) // self._ch
        valid_col = 0 <= col < self._aln.alignment_length
        valid_row = 0 <= row < self._aln.num_sequences
        if valid_col and valid_row:
            self.hovered.emit(row, col)
        if not emit_only:
            self._sel_col = col if valid_col else None
            self._sel_row = row if valid_row else None
            self.update()



# ---------------------------------------------------------------------------
# Public widget
# ---------------------------------------------------------------------------

class AlignmentView(QWidget):
    """Scrollable MSA viewer with color-scheme and zoom controls."""

    def __init__(self, alignment: Optional[SequenceAlignment] = None,
                 parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._aln: Optional[SequenceAlignment] = alignment
        self._setup_ui()
        if alignment:
            self._canvas.set_alignment(alignment)

    def _setup_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Control bar (plain widget — QToolBar clips slider handles)
        bar = QWidget()
        bar.setObjectName("alignment_bar")
        bar.setFixedHeight(44)
        bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        bar_layout = QHBoxLayout(bar)
        bar_layout.setContentsMargins(8, 6, 8, 6)
        bar_layout.setSpacing(6)

        bar_layout.addWidget(QLabel("Color:"))
        self._scheme_cb = QComboBox()
        for name in COLOR_SCHEMES:
            self._scheme_cb.addItem(name)
        self._scheme_cb.currentIndexChanged.connect(self._on_scheme)
        bar_layout.addWidget(self._scheme_cb)

        sep1 = QWidget(); sep1.setFixedWidth(1)
        sep1.setStyleSheet("background: #D1D1D6;")
        bar_layout.addWidget(sep1)

        bar_layout.addWidget(QLabel("Zoom:"))
        self._zoom = QSlider(Qt.Orientation.Horizontal)
        self._zoom.setRange(7, 22)
        self._zoom.setValue(14)
        self._zoom.setFixedWidth(120)
        self._zoom.valueChanged.connect(self._on_zoom)
        bar_layout.addWidget(self._zoom)

        sep2 = QWidget(); sep2.setFixedWidth(1)
        sep2.setStyleSheet("background: #D1D1D6;")
        bar_layout.addWidget(sep2)

        self._info = QLabel("Load an alignment to start")
        self._info.setObjectName("status_label")
        self._info.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        bar_layout.addWidget(self._info, stretch=1)
        layout.addWidget(bar)

        # Scrollable canvas
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(False)
        self._canvas = _AlignmentCanvas()
        self._canvas.hovered.connect(self._on_hover)
        self._scroll.setWidget(self._canvas)
        layout.addWidget(self._scroll)

    def set_alignment(self, aln: SequenceAlignment) -> None:
        self._aln = aln
        self._canvas.set_alignment(aln)

    # ------------------------------------------------------------------

    def _on_scheme(self, _: int) -> None:
        name = self._scheme_cb.currentText()
        self._canvas.set_scheme(COLOR_SCHEMES[name])

    def _on_zoom(self, v: int) -> None:
        self._canvas.set_cell_size(v, int(v * 1.35))

    def _on_hover(self, row: int, col: int) -> None:
        if self._aln:
            rec = self._aln.records[row]
            ch = rec.sequence[col] if col < rec.length else "-"
            self._info.setText(
                f"  {rec.id}  |  pos {col + 1}  |  {ch.upper()}"
            )
