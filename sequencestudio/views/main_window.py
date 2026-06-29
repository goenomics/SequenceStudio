"""Main application window."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QDockWidget, QTabWidget, QLabel,
    QListWidget, QListWidgetItem, QVBoxLayout, QHBoxLayout,
    QFileDialog, QMessageBox, QApplication, QSizePolicy,
    QStatusBar, QFrame, QToolButton, QMenu,
)
from PySide6.QtCore import Qt, QThread, Signal, QSettings, QSize
from PySide6.QtGui import QAction, QKeySequence, QFont, QPixmap

_ICON_PATH = Path(__file__).parent.parent / "assets" / "app_icon.png"

from sequencestudio.models.sequence import SequenceRecord, SequenceAlignment
from sequencestudio.core.io import read_sequences, write_sequences, WRITE_FORMATS
from sequencestudio.core.alignment import align_sequences
from sequencestudio.views.alignment_view import AlignmentView
from sequencestudio.views.blast_dialog import BlastDialog
from sequencestudio.views.restriction_map import RestrictionMapView
from sequencestudio.views.orf_finder import OrfFinderView
from sequencestudio.views.primer_design import PrimerDesignDialog
from sequencestudio.views.dot_plot import DotPlotView
from sequencestudio.views.tree_view import TreeView


# ---------------------------------------------------------------------------
# Welcome screen
# ---------------------------------------------------------------------------

class _WelcomeWidget(QWidget):
    """Styled start-up screen shown in the first tab."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setAutoFillBackground(True)

        outer = QVBoxLayout(self)
        outer.setAlignment(Qt.AlignmentFlag.AlignCenter)
        # Add left border color to match the sidebar
        outer.setContentsMargins(0, 0, 0, 0)


        card = QFrame()
        card.setObjectName("welcome_card")
        card.setStyleSheet(
            "QFrame#welcome_card {"
            "  background: #FFFFFF;"
            "  border: 1px solid #E5E5EA;"
            "  border-radius: 16px;"
            "  padding: 10px;"
            "}"
        )
        card.setFixedSize(420, 360)
        card_layout = QVBoxLayout(card)
        card_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        card_layout.setSpacing(14)

        # Logo
        logo_label = QLabel()
        logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        if _ICON_PATH.exists():
            pix = QPixmap(str(_ICON_PATH)).scaled(
                72, 72,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            logo_label.setPixmap(pix)
        card_layout.addWidget(logo_label)

        title = QLabel("SequenceStudio")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet(
            "font-size: 26px; font-weight: 700; color: #1C1C1E; letter-spacing: -0.5px;"
        )
        card_layout.addWidget(title)

        subtitle = QLabel("Biological sequence analysis workbench")
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setStyleSheet("font-size: 14px; color: #636366;")
        card_layout.addWidget(subtitle)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #E5E5EA;")
        card_layout.addWidget(sep)

        hint = QLabel(
            "Open a sequence file with <b>⌘O</b> to get started.<br/>"
            "<span style='font-size:11px; color:#8E8E93;'>"
            "FASTA · GenBank · Clustal · NEXUS · PHYLIP</span>"
        )
        hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hint.setStyleSheet("font-size: 13px; color: #3C3C43;")
        hint.setTextFormat(Qt.TextFormat.RichText)
        card_layout.addWidget(hint)

        outer.addWidget(card)


# ---------------------------------------------------------------------------
# Background alignment worker
# ---------------------------------------------------------------------------

class _AlignWorker(QThread):
    finished = Signal(object)   # SequenceAlignment
    error    = Signal(str)

    def __init__(self, records: list[SequenceRecord], method: str) -> None:
        super().__init__()
        self._records = records
        self._method  = method

    def run(self) -> None:
        try:
            self.finished.emit(align_sequences(self._records, self._method))
        except Exception as exc:
            self.error.emit(str(exc))


# ---------------------------------------------------------------------------
# Left-panel sequence list
# ---------------------------------------------------------------------------

class _SeqListPanel(QWidget):
    selection_changed = Signal(list)   # list[SequenceRecord]

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._records: list[SequenceRecord] = []
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 6, 4, 4)
        self.setObjectName("seq_panel")

        hdr = QLabel("SEQUENCES")
        hdr.setObjectName("sidebar_header")
        lay.addWidget(hdr)

        self._list = QListWidget()
        self._list.setFont(QFont("Menlo", 10))
        self._list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self._list.setAlternatingRowColors(False)
        self._list.itemSelectionChanged.connect(self._emit_selection)
        lay.addWidget(self._list)

        self._footer = QLabel("")
        self._footer.setObjectName("sidebar_footer")
        lay.addWidget(self._footer)

    def set_records(self, records: list[SequenceRecord]) -> None:
        self._records = records
        self._list.clear()
        for r in records:
            item = QListWidgetItem(f"{r.id}  [{r.length:,} bp]")
            item.setToolTip(
                f"ID: {r.id}\n"
                f"Description: {r.description}\n"
                f"Length: {r.length:,}\n"
                f"Type: {r.molecule_type}\n"
                f"GC: {r.gc_content():.1f}%"
            )
            self._list.addItem(item)
        self._footer.setText(f"{len(records)} sequence(s)")

    def all_records(self) -> list[SequenceRecord]:
        return self._records

    def selected_records(self) -> list[SequenceRecord]:
        rows = [self._list.row(i) for i in self._list.selectedItems()]
        return [self._records[r] for r in rows if r < len(self._records)]

    def _emit_selection(self) -> None:
        self.selection_changed.emit(self.selected_records())


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):

    def __init__(self) -> None:
        super().__init__()
        self._records: list[SequenceRecord] = []
        self._alignment: Optional[SequenceAlignment] = None
        self._worker: Optional[_AlignWorker] = None
        self._settings = QSettings("Goenomics", "SequenceStudio")

        self._build_window()
        self._build_menus()
        self._build_toolbar()
        self._build_panels()
        self._restore_geometry()

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def _build_window(self) -> None:
        self.setWindowTitle("SequenceStudio")
        self.setMinimumSize(1100, 680)
        self.resize(1400, 900)

    def _build_menus(self) -> None:
        mb = self.menuBar()

        # ---- File ----
        fm = mb.addMenu("&File")
        self._a_open = QAction("&Open…", self, shortcut=QKeySequence.StandardKey.Open)
        self._a_open.triggered.connect(self.open_file)
        fm.addAction(self._a_open)

        self._a_save = QAction("Save &Alignment…", self,
                               shortcut=QKeySequence.StandardKey.Save)
        self._a_save.triggered.connect(self.save_alignment)
        self._a_save.setEnabled(False)
        fm.addAction(self._a_save)

        fm.addSeparator()
        a_quit = QAction("&Quit", self, shortcut=QKeySequence.StandardKey.Quit)
        a_quit.triggered.connect(QApplication.quit)
        fm.addAction(a_quit)

        # ---- Edit ----
        em = mb.addMenu("&Edit")
        a_selall = QAction("Select &All", self,
                           shortcut=QKeySequence.StandardKey.SelectAll)
        a_selall.triggered.connect(lambda: self._seq_panel._list.selectAll())
        em.addAction(a_selall)

        # ---- Alignment ----
        am = mb.addMenu("&Alignment")
        self._a_mafft = QAction("Align with &MAFFT", self)
        self._a_mafft.setToolTip("brew install mafft")
        self._a_mafft.triggered.connect(lambda: self._run_alignment("mafft"))
        am.addAction(self._a_mafft)

        self._a_muscle = QAction("Align with M&USCLE", self)
        self._a_muscle.setToolTip("brew install muscle")
        self._a_muscle.triggered.connect(lambda: self._run_alignment("muscle"))
        am.addAction(self._a_muscle)

        self._a_clustalo = QAction("Align with &ClustalOmega", self)
        self._a_clustalo.setToolTip("Source currently unavailable via brew; use MAFFT instead")
        self._a_clustalo.triggered.connect(lambda: self._run_alignment("clustalo"))
        am.addAction(self._a_clustalo)

        am.addSeparator()
        self._a_builtin = QAction("Align (Built-in, no install)", self)
        self._a_builtin.setToolTip("Pure-Python progressive alignment — no external tool needed")
        self._a_builtin.triggered.connect(lambda: self._run_alignment("builtin"))
        am.addAction(self._a_builtin)

        # ---- Analysis ----
        anlm = mb.addMenu("A&nalysis")

        self._a_blast = QAction("&BLAST Search…", self)
        self._a_blast.triggered.connect(self._open_blast)
        anlm.addAction(self._a_blast)

        self._a_re = QAction("&Restriction Map…", self)
        self._a_re.triggered.connect(self._open_restriction)
        anlm.addAction(self._a_re)

        self._a_orf = QAction("&ORF Finder…", self)
        self._a_orf.triggered.connect(self._open_orf)
        anlm.addAction(self._a_orf)

        self._a_primer = QAction("&Primer Design…", self)
        self._a_primer.triggered.connect(self._open_primer)
        anlm.addAction(self._a_primer)

        anlm.addSeparator()

        self._a_dot = QAction("&Dot Plot…", self)
        self._a_dot.triggered.connect(self._open_dotplot)
        anlm.addAction(self._a_dot)

        self._a_tree = QAction("Phylogenetic &Tree…", self)
        self._a_tree.triggered.connect(self._open_tree)
        anlm.addAction(self._a_tree)

    def _build_toolbar(self) -> None:
        tb = self.addToolBar("Main")
        tb.setMovable(False)
        tb.setIconSize(QSize(16, 16))

        tb.addAction(self._a_open)
        tb.addSeparator()

        # ---- Alignment split-button dropdown ----
        self._align_method = "mafft"
        self._align_btn = QToolButton()
        self._align_btn.setText("Align: MAFFT")
        self._align_btn.setToolTip("Run multiple sequence alignment")
        self._align_btn.setPopupMode(
            QToolButton.ToolButtonPopupMode.MenuButtonPopup
        )
        self._align_btn.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextOnly
        )
        self._align_btn.setMinimumWidth(130)
        self._align_btn.clicked.connect(
            lambda: self._run_alignment(self._align_method)
        )

        _ALIGN_METHODS = [
            ("MAFFT",    "mafft",   "Recommended  ·  brew install mafft"),
            ("MUSCLE",   "muscle",  "brew install muscle"),
            ("ClustalΩ", "clustalo", "Source currently unavailable via brew; use MAFFT instead"),
            ("Built-in", "builtin", "Pure-Python fallback — no external tool required"),
        ]
        align_menu = QMenu(self._align_btn)
        for label, key, tip in _ALIGN_METHODS:
            act = align_menu.addAction(label)
            act.setToolTip(tip)
            act.triggered.connect(
                lambda checked=False, k=key, l=label: self._set_align_method(k, l)
            )
        self._align_btn.setMenu(align_menu)
        tb.addWidget(self._align_btn)
        tb.addSeparator()

        for a in (
            self._a_blast, self._a_re, self._a_orf, self._a_primer, None,
            self._a_dot, self._a_tree,
        ):
            if a is None:
                tb.addSeparator()
            else:
                tb.addAction(a)

    def _set_align_method(self, method: str, label: str) -> None:
        """Update the active alignment method shown on the toolbar button."""
        self._align_method = method
        self._align_btn.setText(f"Align: {label}")

    def _build_panels(self) -> None:
        # Left dock
        self._seq_panel = _SeqListPanel()
        self._seq_panel.selection_changed.connect(self._on_selection)


        dock = QDockWidget("Sequences", self)
        dock.setWidget(self._seq_panel)
        dock.setMinimumWidth(210)
        dock.setAllowedAreas(
            Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea
        )
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

        # Central tabs
        self._tabs = QTabWidget()
        self._tabs.setTabsClosable(True)
        self._tabs.tabCloseRequested.connect(self._close_tab)
        self._tabs.setDocumentMode(True)
        self._tabs.tabBar().setDrawBase(False)

        welcome = _WelcomeWidget()
        self._tabs.addTab(welcome, "Start")
        # make welcome tab non-closable
        self._tabs.tabBar().setTabButton(
            0, self._tabs.tabBar().ButtonPosition.RightSide, None
        )

        self.setCentralWidget(self._tabs)

    def _restore_geometry(self) -> None:
        geo = self._settings.value("geometry")
        if geo:
            self.restoreGeometry(geo)

    def closeEvent(self, event) -> None:
        self._settings.setValue("geometry", self.saveGeometry())
        super().closeEvent(event)

    # ------------------------------------------------------------------
    # File operations
    # ------------------------------------------------------------------

    def open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "",
            "All Supported (*.fasta *.fa *.fna *.faa *.gb *.gbk *.aln *.nex *.phy);;"
            "FASTA (*.fasta *.fa *.fna *.faa);;"
            "GenBank (*.gb *.gbk);;"
            "Clustal (*.aln);;"
            "NEXUS (*.nex);;"
            "PHYLIP (*.phy);;"
            "All Files (*)",
        )
        if not path:
            return
        try:
            records = read_sequences(path)
        except Exception as exc:
            QMessageBox.critical(self, "Open Error", str(exc))
            return

        self._records = records
        self._seq_panel.set_records(records)
        fname = Path(path).name
        self.setWindowTitle(f"SequenceStudio — {fname}")
        self.statusBar().showMessage(f"Loaded {len(records)} sequence(s) from {fname}")

        # Auto-show alignment if all sequences have the same length
        if len(records) > 1 and len({r.length for r in records}) == 1:
            aln = SequenceAlignment(records=records, name=Path(path).stem)
            self._alignment = aln
            self._show_alignment(aln)
            self._a_save.setEnabled(True)

    def save_alignment(self) -> None:
        if not self._alignment:
            return
        path, sel = QFileDialog.getSaveFileName(
            self, "Save Alignment", "",
            "FASTA (*.fasta);;Clustal (*.aln);;NEXUS (*.nex);;PHYLIP (*.phy)",
        )
        if not path:
            return
        fmt_map = {"FASTA": "fasta", "Clustal": "clustal", "NEXUS": "nexus", "PHYLIP": "phylip"}
        fmt = fmt_map.get(sel.split(" ")[0], "fasta")
        try:
            write_sequences(path, self._alignment.records, fmt)
            self.statusBar().showMessage(f"Saved alignment to {Path(path).name}")
        except Exception as exc:
            QMessageBox.critical(self, "Save Error", str(exc))

    # ------------------------------------------------------------------
    # Tab management
    # ------------------------------------------------------------------

    def _close_tab(self, idx: int) -> None:
        if self._tabs.count() > 1:
            self._tabs.removeTab(idx)

    def _add_tab(self, widget: QWidget, title: str) -> None:
        idx = self._tabs.addTab(widget, title)
        self._tabs.setCurrentIndex(idx)

    def _show_alignment(self, aln: SequenceAlignment) -> None:
        view = AlignmentView(aln)
        self._add_tab(view, f"Aln: {aln.name}")

    # ------------------------------------------------------------------
    # Status updates
    # ------------------------------------------------------------------

    def _on_selection(self, records: list[SequenceRecord]) -> None:
        if not records:
            self.statusBar().showMessage("Ready")
        elif len(records) == 1:
            r = records[0]
            self.statusBar().showMessage(
                f"{r.id}  |  {r.length:,} bp  |  {r.molecule_type}  |  GC {r.gc_content():.1f}%"
            )
        else:
            self.statusBar().showMessage(f"{len(records)} sequences selected")

    # ------------------------------------------------------------------
    # Alignment
    # ------------------------------------------------------------------

    def _run_alignment(self, method: str) -> None:
        records = self._seq_panel.selected_records() or self._seq_panel.all_records()
        if len(records) < 2:
            QMessageBox.information(self, "Alignment", "Select at least 2 sequences.")
            return
        self.setEnabled(False)
        self.statusBar().showMessage(f"Running {method} alignment on {len(records)} sequences…")
        self._worker = _AlignWorker(records, method)
        self._worker.finished.connect(self._on_align_done)
        self._worker.error.connect(self._on_align_error)
        self._worker.start()

    def _on_align_done(self, aln: SequenceAlignment) -> None:
        self.setEnabled(True)
        self._alignment = aln
        self._seq_panel.set_records(aln.records)
        self._show_alignment(aln)
        self._a_save.setEnabled(True)
        self.statusBar().showMessage(
            f"Alignment done: {aln.num_sequences} seqs × {aln.alignment_length} positions"
        )

    def _on_align_error(self, msg: str) -> None:
        self.setEnabled(True)
        QMessageBox.critical(self, "Alignment Error", msg)
        self.statusBar().showMessage("Alignment failed")

    # ------------------------------------------------------------------
    # Analysis actions
    # ------------------------------------------------------------------

    def _active_record(self) -> Optional[SequenceRecord]:
        sel = self._seq_panel.selected_records()
        return sel[0] if sel else (self._seq_panel.all_records() or [None])[0]

    def _require_record(self, tool_name: str) -> Optional[SequenceRecord]:
        rec = self._active_record()
        if not rec:
            QMessageBox.information(self, tool_name, "Load a sequence file first.")
        return rec

    def _open_blast(self) -> None:
        rec = self._require_record("BLAST")
        if rec:
            BlastDialog(rec, self).show()

    def _open_restriction(self) -> None:
        rec = self._require_record("Restriction Map")
        if rec:
            self._add_tab(RestrictionMapView(rec), f"RE: {rec.id}")

    def _open_orf(self) -> None:
        rec = self._require_record("ORF Finder")
        if rec:
            self._add_tab(OrfFinderView(rec), f"ORFs: {rec.id}")

    def _open_primer(self) -> None:
        rec = self._require_record("Primer Design")
        if rec:
            PrimerDesignDialog(rec, self).show()

    def _open_dotplot(self) -> None:
        records = self._seq_panel.all_records()
        if not records:
            QMessageBox.information(self, "Dot Plot", "Load at least one sequence file first.")
            return
        self._add_tab(DotPlotView(records), "Dot Plot")

    def _open_tree(self) -> None:
        if not self._alignment or self._alignment.num_sequences < 3:
            QMessageBox.information(
                self, "Phylogenetic Tree",
                "An alignment with ≥ 3 sequences is required.\n"
                "Run an alignment first.",
            )
            return
        self._add_tab(TreeView(self._alignment), "Tree")
