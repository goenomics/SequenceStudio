"""Global application theme — macOS-inspired light design on Qt Fusion base."""

from __future__ import annotations

from PySide6.QtGui import QPalette, QColor, QFont
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt

# ---------------------------------------------------------------------------
# Colour tokens
# ---------------------------------------------------------------------------
_C = {
    "bg":           "#FFFFFF",
    "bg_panel":     "#F5F5F7",
    "bg_hover":     "#E8E8ED",
    "bg_press":     "#D1D1D6",
    "bg_input":     "#FFFFFF",
    "accent":       "#007AFF",
    "accent_hover": "#0A84FF",
    "accent_press": "#0060DF",
    "green":        "#34C759",
    "green_hover":  "#2DB84E",
    "border":       "#D1D1D6",
    "border_light": "#E5E5EA",
    "text":         "#1C1C1E",
    "text_sec":     "#636366",
    "text_dis":     "#AEAEB2",
    "scrollbar":    "#C7C7CC",
}

# ---------------------------------------------------------------------------
# Full application stylesheet
# ---------------------------------------------------------------------------
STYLESHEET = f"""
/* ════════════════════════════════════════════════════════
   SequenceStudio · Modern Light Theme
   ════════════════════════════════════════════════════════ */

/* Base ─────────────────────────────────────────────────── */
QMainWindow, QDialog, QWidget {{
    background-color: {_C['bg']};
    color: {_C['text']};
    font-size: 13px;
}}

/* Sequence panel right border ────────────────────── */
QWidget#seq_panel {{
    border-right: 1px solid {_C['border']};
}}

/* Dock resize handle ────────────────────────────────────── */
QMainWindow::separator {{
    background: {_C['bg_panel']};
    width: 4px;
    height: 4px;
}}
QMainWindow::separator:hover {{
    background: {_C['accent']};
}}

/* Menu bar ─────────────────────────────────────────────── */
QMenuBar {{
    background-color: {_C['bg_panel']};
    border-bottom: 1px solid {_C['border']};
    padding: 1px 4px;
    spacing: 0;
}}
QMenuBar::item {{
    padding: 4px 10px;
    border-radius: 5px;
}}
QMenuBar::item:selected, QMenuBar::item:pressed {{
    background: {_C['accent']};
    color: #FFFFFF;
}}

QMenu {{
    background: {_C['bg']};
    border: 1px solid {_C['border']};
    border-radius: 10px;
    padding: 5px 0;
}}
QMenu::item {{
    padding: 5px 18px 5px 14px;
    border-radius: 5px;
    margin: 1px 4px;
}}
QMenu::item:selected {{
    background: {_C['accent']};
    color: #FFFFFF;
}}
QMenu::separator {{
    height: 1px;
    background: {_C['border_light']};
    margin: 4px 12px;
}}

/* Tool bar ─────────────────────────────────────────────── */
QToolBar {{
    background: {_C['bg_panel']};
    border: none;
    border-bottom: 1px solid {_C['border']};
    spacing: 2px;
    padding: 3px 6px;
}}
QToolBar::separator {{
    background: {_C['border']};
    width: 1px;
    margin: 3px 5px;
}}
QToolButton {{
    background: transparent;
    border: none;
    border-radius: 6px;
    padding: 4px 10px;
    color: {_C['text']};
    font-size: 12px;
}}
QToolButton:hover  {{ background: {_C['bg_hover']};  }}
QToolButton:pressed {{ background: {_C['bg_press']};  }}

/* Alignment view control bar ───────────────────────────── */
QWidget#alignment_bar {{
    background: {_C['bg']};
    border-bottom: 1px solid {_C['border']};
}}
QWidget#alignment_bar QLabel {{
    font-size: 12px;
    color: {_C['text_sec']};
}}

/* Dock widget ──────────────────────────────────────────── */
QDockWidget {{
    /* border: none; */
    border-right: 1px solid {_C['border']};
}}
QDockWidget::title {{
    background: {_C['bg_panel']};
    text-align: left;
    padding: 6px 10px;
    border-bottom: 1px solid {_C['border']};
    border-right: 1px solid {_C['border']};
    font-size: 11px;
    font-weight: 700;
    color: {_C['text_sec']};
    letter-spacing: 0.5px;
    text-transform: uppercase;
}}

/* Tab widget ────────────────────────────────────────────── */
QTabWidget {{
    border: none;
}}
QTabWidget::pane {{
    background: {_C['bg']};
    border: none;
}}
QTabBar {{
    background: {_C['bg']};
    border: none;
}}
QTabBar::tab {{
    background: {_C['bg_panel']};
    padding: 7px 20px;
    font-size: 12px;
    color: {_C['text_sec']};
    border: 1px solid {_C['border']};
    border-top-left-radius: 7px;
    border-top-right-radius: 7px;
    border-bottom-left-radius: 0px;
    border-bottom-right-radius: 0px;
    margin-top: 2px;
}}
QTabBar::tab:selected {{
    background: {_C['bg']};
    color: {_C['accent']};
    border-bottom: 1px solid {_C['bg']};
    font-weight: 600;
}}
QTabBar::tab:hover:!selected {{
    background: {_C['bg_hover']};
    color: {_C['text']};
}}

/* Push buttons ─────────────────────────────────────────── */
QPushButton {{
    background: {_C['bg_panel']};
    border: 1px solid {_C['border']};
    border-radius: 7px;
    padding: 5px 16px;
    font-size: 13px;
    color: {_C['text']};
    min-height: 22px;
}}
QPushButton:hover   {{ background: {_C['bg_hover']};  border-color: #AEAEB2; }}
QPushButton:pressed {{ background: {_C['bg_press']};  }}
QPushButton:disabled {{
    color: {_C['text_dis']};
    border-color: {_C['border_light']};
    background: #FAFAFA;
}}

/* Combo box ─────────────────────────────────────────────── */
QComboBox {{
    background: {_C['bg_input']};
    border: 1px solid {_C['border']};
    border-radius: 7px;
    padding: 4px 10px;
    font-size: 13px;
    min-height: 22px;
    selection-background-color: {_C['accent']};
}}
QComboBox:hover {{ border-color: #AEAEB2; }}
QComboBox:focus {{ border-color: {_C['accent']}; }}
QComboBox::drop-down {{ border: none; width: 20px; }}
QComboBox QAbstractItemView {{
    background: {_C['bg']};
    border: 1px solid {_C['border']};
    border-radius: 8px;
    selection-background-color: {_C['accent']};
    selection-color: #FFFFFF;
    padding: 3px;
    outline: none;
}}

/* Spin boxes ────────────────────────────────────────────── */
QSpinBox, QDoubleSpinBox {{
    background: {_C['bg_input']};
    border: 1px solid {_C['border']};
    border-radius: 7px;
    padding: 4px 8px;
    font-size: 13px;
    selection-background-color: {_C['accent']};
}}
QSpinBox:focus, QDoubleSpinBox:focus {{ border-color: {_C['accent']}; }}
QSpinBox::up-button, QDoubleSpinBox::up-button,
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    border: none;
    background: transparent;
    width: 16px;
}}

/* Check box ─────────────────────────────────────────────── */
QCheckBox {{ spacing: 7px; font-size: 13px; }}
QCheckBox::indicator {{
    width: 16px;
    height: 16px;
    border: 1.5px solid {_C['text_dis']};
    border-radius: 4px;
    background: {_C['bg']};
}}
QCheckBox::indicator:hover   {{ border-color: {_C['text_sec']}; }}
QCheckBox::indicator:checked {{
    background: {_C['accent']};
    border-color: {_C['accent']};
}}

/* Group box ─────────────────────────────────────────────── */
QGroupBox {{
    font-size: 11px;
    font-weight: 700;
    color: {_C['text_sec']};
    letter-spacing: 0.4px;
    text-transform: uppercase;
    border: 1px solid {_C['border_light']};
    border-radius: 10px;
    margin-top: 10px;
    padding: 18px 8px 8px 8px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    top: 0px;
    padding: 0 4px;
    background: {_C['bg']};
}}

/* List widget ───────────────────────────────────────────── */
QListWidget {{
    background: transparent;
    border: none;
    outline: none;
    font-size: 12px;
}}
QListWidget::item {{
    padding: 5px 10px;
    border-radius: 6px;
    margin: 1px 4px;
    color: {_C['text']};
}}
QListWidget::item:selected {{
    background: {_C['accent']};
    color: #FFFFFF;
}}
QListWidget::item:hover:!selected {{
    background: {_C['bg_hover']};
}}

/* Table widget ──────────────────────────────────────────── */
QTableWidget {{
    background: {_C['bg']};
    border: none;
    gridline-color: {_C['bg_panel']};
    outline: none;
    font-size: 12px;
}}
QTableWidget::item {{
    padding: 4px 8px;
    border: none;
}}
QTableWidget::item:selected {{
    background: {_C['accent']};
    color: #FFFFFF;
}}
QHeaderView::section {{
    background: {_C['bg_panel']};
    color: {_C['text_sec']};
    font-size: 11px;
    font-weight: 700;
    padding: 5px 8px;
    border: none;
    border-bottom: 1px solid {_C['border']};
    border-right: 1px solid {_C['border_light']};
    letter-spacing: 0.3px;
    text-transform: uppercase;
}}
QHeaderView::section:last {{ border-right: none; }}

/* Scroll bars ───────────────────────────────────────────── */
QScrollBar:vertical {{
    background: transparent;
    width: 8px;
    margin: 0;
    border: none;
}}
QScrollBar::handle:vertical {{
    background: {_C['scrollbar']};
    border-radius: 4px;
    min-height: 28px;
}}
QScrollBar::handle:vertical:hover {{ background: #AEAEB2; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; border: none; }}
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {{ background: transparent; }}

QScrollBar:horizontal {{
    background: transparent;
    height: 8px;
    margin: 0;
    border: none;
}}
QScrollBar::handle:horizontal {{
    background: {_C['scrollbar']};
    border-radius: 4px;
    min-width: 28px;
}}
QScrollBar::handle:horizontal:hover {{ background: #AEAEB2; }}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{ width: 0; border: none; }}
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {{ background: transparent; }}

/* Text edit ─────────────────────────────────────────────── */
QTextEdit, QPlainTextEdit {{
    background: {_C['bg']};
    border: 1px solid {_C['border_light']};
    border-radius: 8px;
    padding: 6px;
    selection-background-color: {_C['accent']};
    color: {_C['text']};
}}
QTextEdit:focus, QPlainTextEdit:focus {{ border-color: {_C['accent']}; }}

/* Status bar ────────────────────────────────────────────── */
QStatusBar {{
    background: {_C['bg_panel']};
    border-top: 1px solid {_C['border']};
    color: {_C['text_sec']};
    font-size: 11px;
    padding: 2px 8px;
}}
QStatusBar::item {{ border: none; }}

/* Splitter ──────────────────────────────────────────────── */
QSplitter::handle {{ background: {_C['border']}; }}
QSplitter::handle:horizontal {{ width: 1px; }}
QSplitter::handle:vertical   {{ height: 1px; }}

/* Scroll area ───────────────────────────────────────────── */
QScrollArea {{ border: none; background: transparent; }}

/* Slider ────────────────────────────────────────────────── */
QSlider::groove:horizontal {{
    background: {_C['border_light']};
    border-radius: 2px;
    height: 3px;
}}
QSlider::sub-page:horizontal {{
    background: {_C['accent']};
    border-radius: 2px;
}}
QSlider::handle:horizontal {{
    background: #FFFFFF;
    border: 2px solid {_C['accent']};
    border-radius: 5px;
    width: 10px;
    height: 10px;
    margin: -4px 0;
}}

/* Progress bar ──────────────────────────────────────────── */
QProgressBar {{
    background: {_C['border_light']};
    border: none;
    border-radius: 4px;
    height: 5px;
    text-align: center;
    font-size: 10px;
    color: {_C['text_sec']};
}}
QProgressBar::chunk {{
    background: {_C['accent']};
    border-radius: 4px;
}}

/* Tool tips ─────────────────────────────────────────────── */
QToolTip {{
    background: {_C['text']};
    color: #FFFFFF;
    border: none;
    border-radius: 6px;
    padding: 5px 9px;
    font-size: 12px;
    opacity: 230;
}}

/* ── Semantic object names ───────────────────────────────── */
QLabel#status_label {{
    color: {_C['text_sec']};
    font-style: italic;
    font-size: 12px;
}}
QLabel#sidebar_header {{
    color: {_C['text_sec']};
    font-size: 10px;
    font-weight: 700;
    letter-spacing: 0.8px;
    text-transform: uppercase;
    padding: 0 2px;
}}
QLabel#sidebar_footer {{
    color: {_C['text_dis']};
    font-size: 10px;
    padding: 0 4px;
}}
"""


def build_palette() -> QPalette:
    """Return a QPalette that matches the light theme."""
    pal = QPalette()
    c = QColor

    pal.setColor(QPalette.ColorRole.Window,          c(_C["bg_panel"]))
    pal.setColor(QPalette.ColorRole.WindowText,       c(_C["text"]))
    pal.setColor(QPalette.ColorRole.Base,             c(_C["bg"]))
    pal.setColor(QPalette.ColorRole.AlternateBase,    c("#F9F9F9"))
    pal.setColor(QPalette.ColorRole.ToolTipBase,      c(_C["text"]))
    pal.setColor(QPalette.ColorRole.ToolTipText,      c("#FFFFFF"))
    pal.setColor(QPalette.ColorRole.Text,             c(_C["text"]))
    pal.setColor(QPalette.ColorRole.Button,           c(_C["bg_panel"]))
    pal.setColor(QPalette.ColorRole.ButtonText,       c(_C["text"]))
    pal.setColor(QPalette.ColorRole.BrightText,       c("#FFFFFF"))
    pal.setColor(QPalette.ColorRole.Link,             c(_C["accent"]))
    pal.setColor(QPalette.ColorRole.Highlight,        c(_C["accent"]))
    pal.setColor(QPalette.ColorRole.HighlightedText,  c("#FFFFFF"))
    pal.setColor(QPalette.ColorRole.Midlight,         c(_C["bg_hover"]))
    pal.setColor(QPalette.ColorRole.Mid,              c(_C["border"]))
    pal.setColor(QPalette.ColorRole.Dark,             c(_C["text_dis"]))
    pal.setColor(QPalette.ColorRole.Shadow,           c(_C["text_sec"]))

    # Disabled state
    pal.setColor(QPalette.ColorGroup.Disabled,
                 QPalette.ColorRole.WindowText, c(_C["text_dis"]))
    pal.setColor(QPalette.ColorGroup.Disabled,
                 QPalette.ColorRole.Text,        c(_C["text_dis"]))
    pal.setColor(QPalette.ColorGroup.Disabled,
                 QPalette.ColorRole.ButtonText,  c(_C["text_dis"]))

    return pal


def apply_theme(app: QApplication) -> None:
    """Apply Fusion + custom palette + full QSS to *app*."""
    app.setStyle("Fusion")
    app.setPalette(build_palette())

    # System font with explicit fallback chain
    font = QFont("-apple-system")
    if not font.exactMatch():
        font = QFont("SF Pro Text")
    if not font.exactMatch():
        font = QFont("Helvetica Neue")
    font.setPointSize(13)
    font.setStyleHint(QFont.StyleHint.SansSerif)
    app.setFont(font)

    app.setStyleSheet(STYLESHEET)
