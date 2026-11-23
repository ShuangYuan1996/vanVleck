"""vanVleck entry point: python main.py [file1.cif file2.xyz ...]"""
import sys

from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from sco.gui.main_window import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("vanVleck")
    app.setFont(QFont("Segoe UI", 9))
    win = MainWindow()
    win.show()
    if len(sys.argv) > 1:
        win.load_many(sys.argv[1:])     # several files may be given
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
