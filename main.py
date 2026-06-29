"""SequenceStudio — Entry point."""

import sys
import os

# macOS layer-backing required for Qt6
os.environ.setdefault("QT_MAC_WANTS_LAYER", "1")

from sequencestudio.app import SequenceStudioApp


def main():
    app = SequenceStudioApp(sys.argv)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
