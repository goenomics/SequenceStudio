"""SequenceStudio application bootstrap."""

import sys
from pathlib import Path
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon

from sequencestudio.utils.style import apply_theme

_ICON_PATH = Path(__file__).parent / "assets" / "app_icon.png"


class SequenceStudioApp(QApplication):
    def __init__(self, argv: list[str]):
        super().__init__(argv)
        self.setApplicationName("SequenceStudio")
        self.setApplicationVersion("1.0.0")
        self.setOrganizationName("Goenomics")
        self.setOrganizationDomain("goenomics.com")

        if _ICON_PATH.exists():
            self.setWindowIcon(QIcon(str(_ICON_PATH)))

        apply_theme(self)

        # Defer main window import to avoid circular imports at module level
        from sequencestudio.views.main_window import MainWindow
        self._main_window = MainWindow()
        self._main_window.show()
