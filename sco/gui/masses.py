"""Input widget for the masses used in the T1u coordinates (metal M0, ligand atom M)."""
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QWidget

from .. import vanvleck as vv


class MassEditor(QWidget):
    changed = Signal(float, float)        # (M0 metal, M ligand atom)

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        tip = ("Masses used in the T1u coordinates Q10–Q15 (only the ratio M0/M matters).\n"
               "The six donor atoms are treated as identical atoms of mass M.")
        self.e0 = QLineEdit(f"{vv.DEFAULT_M0:g}")
        self.e1 = QLineEdit(f"{vv.DEFAULT_M:g}")
        for lab, e in (("Metal mass M₀:", self.e0), ("Ligand mass M:", self.e1)):
            l_ = QLabel(lab)
            l_.setToolTip(tip)
            e.setToolTip(tip)
            e.setValidator(QDoubleValidator(1e-6, 1e6, 6))
            e.setMaximumWidth(70)
            e.editingFinished.connect(self._emit)
            lay.addWidget(l_)
            lay.addWidget(e)
        lay.addStretch(1)
        self._last = (vv.DEFAULT_M0, vv.DEFAULT_M)

    def values(self) -> tuple[float, float]:
        out = []
        for e, d in ((self.e0, vv.DEFAULT_M0), (self.e1, vv.DEFAULT_M)):
            try:
                v = float(e.text())
                out.append(v if v > 0 else d)
            except ValueError:
                out.append(d)
        return out[0], out[1]

    def set_values(self, m0: float, m: float):
        self.e0.setText(f"{m0:g}")
        self.e1.setText(f"{m:g}")
        self._last = (m0, m)

    def _emit(self):
        v = self.values()
        if v != self._last:
            self._last = v
            self.changed.emit(*v)
