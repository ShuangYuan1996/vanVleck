"""Compare window: molecules from any structures are collected in one list.

- Same material: the molecule with the highest inner-sphere symmetry is the reference (first added
  wins a tie); the frames of all molecules are chosen jointly among the equivalent standard frames
  so that Q7 has one sign and the Q values are as close as possible, making Q directly comparable.
- Different chemical structures: a warning is shown and each Q uses its own best frame.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QGuiApplication
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QHBoxLayout, QHeaderView, QLabel, QPushButton,
                               QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .. import symmetry as sym
from .. import vanvleck as vv
from .. import distortion as dist
from ..coordination import CoordinationSphere
from ..materials import formula_text, molecule_key, same_material
from .masses import MassEditor
from .qchart import SERIES_COLORS, QBarChart, export_q_dialog, fmt3


@dataclass
class CompareEntry:
    name: str                     # legend name
    dataset: str
    temperature: float | None
    metal: str
    coords: np.ndarray            # (7,3)
    elems: list
    key: tuple                    # chemical identity of the molecule
    symres: sym.SymmetryResult
    result: vv.VanVleckResult | None = None
    role: str = ""
    dist: dist.Distortion | None = None      # zeta, Delta, Sigma, Theta, V (display only; not exported)


class CompareDialog(QDialog):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Compare van Vleck coordinates")
        self.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        self.resize(1100, 760)
        self.entries: list[CompareEntry] = []
        self.tol = sym.DEFAULT_TOL
        self.m0, self.m = vv.DEFAULT_M0, vv.DEFAULT_M

        lay = QVBoxLayout(self)
        self.msg = QLabel("Click complexes in the main window while this window is open; each one is added here.")
        self.msg.setWordWrap(True)
        self.msg.setTextFormat(Qt.RichText)
        lay.addWidget(self.msg)
        self.mass_edit = MassEditor()
        lay.addWidget(self.mass_edit)

        split = QSplitter(Qt.Vertical)
        top = QWidget()
        tl = QVBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(0, 15)
        self.table.setHorizontalHeaderLabels(["#", "Dataset", "T (K)", "Metal", "Symmetry", "|Q| (Å)",
                                              "ζ (Å)", "Δ", "Σ (°)", "Θ (°)", "Θmin (°)", "Θmax (°)",
                                              "V (Å³)", "RMSD (Å)", "Frame"])
        self.table.horizontalHeaderItem(13).setToolTip(
            "RMSD of the 7-atom inner sphere from the reference after best superposition (Kabsch);\n"
            "atoms are matched through the aligned standard frames. Same material only.")
        self.ref_entry = None
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        tl.addWidget(self.table)
        row = QHBoxLayout()
        b_rm = QPushButton("Remove selected")
        b_rm.clicked.connect(self.remove_selected)
        b_clear = QPushButton("Clear all")
        b_clear.clicked.connect(self.clear)
        b_copy = QPushButton("Copy Q table")
        b_copy.clicked.connect(self.copy_table)
        b_exp = QPushButton("Export .txt")
        b_exp.setToolTip("Save Q1–Q15 of all molecules as a tab-separated .txt file")
        b_exp.clicked.connect(self.export_txt)
        b_img = QPushButton("Copy image")
        b_img.clicked.connect(lambda: QGuiApplication.clipboard().setPixmap(self.chart.grab()))
        for b in (b_rm, b_clear):
            row.addWidget(b)
        row.addStretch(1)
        for b in (b_exp, b_copy, b_img):
            row.addWidget(b)
        tl.addLayout(row)
        split.addWidget(top)
        self.chart = QBarChart()
        split.addWidget(self.chart)
        split.setSizes([200, 560])
        lay.addWidget(split, 1)

    # ------------------------------------------------------------------ data
    def add(self, sphere: CoordinationSphere, symres: sym.SymmetryResult, dataset: str,
            temperature: float | None, tol: float):
        """Add a molecule; return (entry, newly_added). Duplicates are not added again."""
        if tol != self.tol:
            self.tol = tol
            if self.isVisible():
                self.recompute()
        for e in self.entries:
            if e.dataset == dataset and e.metal == sphere.metal.label and np.allclose(e.coords, sphere.coordinates()):
                return e, False
        T = f" {temperature:g} K" if temperature is not None else ""
        elems = [sphere.metal.element] + [a.element for a in sphere.ligands]
        self.entries.append(CompareEntry(f"{dataset}{T} {sphere.metal.label}", dataset, temperature,
                                         sphere.metal.label, sphere.coordinates(), elems,
                                         molecule_key(sphere), symres,
                                         dist=dist.compute(sphere.coordinates())))
        if self.isVisible():           # compute only while the window is visible
            self.recompute()
        return self.entries[-1], True

    def remove_selected(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            del self.entries[r]
        self.recompute()

    def clear(self):
        self.entries.clear()
        self.recompute()

    def set_masses(self, m0: float, m: float):
        self.m0, self.m = m0, m
        self.mass_edit.set_values(m0, m)
        if self.isVisible():
            self.recompute()

    def _cands(self, e):
        return vv.compute(e.coords, e.elems, e.symres, zero_tol=self.tol,
                          m_metal=self.m0, m_lig=self.m, candidates_only=True)

    def _joint_reference_frame(self, ref, others) -> int:
        """Among the equivalent frames of the reference, pick the one for which the other molecules can
        all have the same sign of Q7, then the closest Q; remaining ties prefer Q7 > 0.
        Returns its index (the `pick` argument of vv.compute)."""
        ref_c = self._cands(ref)
        oth_c = [self._cands(e) for e in others]

        def sgn(q7):
            return 0 if abs(q7) < self.tol else (1 if q7 > 0 else -1)

        best, best_key = 0, None
        for k, Qr in enumerate(ref_c):
            s_ref = sgn(Qr[6])
            mism, dsum = 0, 0.0
            for cands in oth_c:
                opts = [(1 if (s_ref and sgn(Q[6]) and sgn(Q[6]) != s_ref) else 0,
                         float(np.linalg.norm(Q[1:] - Qr[1:]))) for Q in cands]
                m_, d_ = min(opts)
                mism += m_
                dsum += d_
            key = (mism, round(dsum, 9), 0 if s_ref >= 0 else 1, k)
            if best_key is None or key < best_key:
                best, best_key = k, key
        return best

    def _q7_sign(self, result) -> int:
        q7 = result.Q[6]
        return 0 if abs(q7) < self.tol else int(np.sign(q7))

    def recompute(self):
        E = self.entries
        if not E:
            self.msg.setText("Click complexes in the main window while this window is open; each one is added here.")
            self._fill_table()
            self.chart.set_series([], [], vv.IRREPS, self.tol, "")
            self.changed.emit()
            return
        same = all(same_material(E[0].key, e.key) for e in E[1:])
        self.ref_entry = None
        for e in E:
            e.result = None
        ok_entries = [e for e in E if e.symres.label in vv.STANDARD_SUBGROUPS]
        if same and ok_entries:
            ref = max(ok_entries, key=lambda e: e.symres.order)      # on a tie, max returns the first added
            k = self._joint_reference_frame(ref, [e for e in ok_entries if e is not ref])
            ref.result = vv.compute(ref.coords, ref.elems, ref.symres, zero_tol=self.tol,
                                    m_metal=self.m0, m_lig=self.m, pick=k)
            ref.role = "reference"
            self.ref_entry = ref
            for e in E:
                if e is ref:
                    continue
                if e.symres.label not in vv.STANDARD_SUBGROUPS:
                    e.result, e.role = None, "no standard frame"
                    continue
                e.result = vv.compute(e.coords, e.elems, e.symres, zero_tol=self.tol, ref_Q=ref.result.Q,
                                      m_metal=self.m0, m_lig=self.m, q7_sign=self._q7_sign(ref.result))
                e.role = f"aligned to #{E.index(ref) + 1}"
            self.msg.setText(f"<span style='color:#080'>✔ All {len(E)} molecules are the same material "
                             f"({formula_text(E[0].key)}).</span> Reference: <b>#{E.index(ref) + 1} "
                             f"{ref.name}</b> ({sym.html(ref.symres.label)}, highest symmetry); the other "
                             f"frames are chosen among the equivalent standard frames to match it.")
        else:
            for e in E:                  # different materials: own best frame, Q7 > 0 where possible
                if e.symres.label in vv.STANDARD_SUBGROUPS:
                    e.result = vv.compute(e.coords, e.elems, e.symres, zero_tol=self.tol,
                                          m_metal=self.m0, m_lig=self.m, q7_sign=1)
                    e.role = "own best frame"
                else:
                    e.result, e.role = None, "no standard frame"
            kinds = sorted({formula_text(e.key) for e in E})
            self.msg.setText("<span style='color:#c00'><b>✘ These molecules are not the same material and "
                             "cannot be compared.</b></span> Each Q is calculated in its own best frame.<br>"
                             + "<br>".join(f"&nbsp;&nbsp;• {k}" for k in kinds))
        self._fill_table()
        shown = [e for e in E if e.result is not None]
        title = "van Vleck coordinates" + ("" if same else " (different materials — independent frames)")
        self.chart.set_series([e.result.Q for e in shown],
                              [f"#{E.index(e) + 1} {e.name}\n|Q| = {e.result.norm:.4f} Å   ζ = {e.dist.zeta:.4f} Å   "
                               f"Δ = {dist.sci(e.dist.delta)}   Σ = {e.dist.sigma:.1f}°   Θ = {e.dist.theta:.1f}°"
                               for e in shown],
                              vv.IRREPS, self.tol, title)
        self.changed.emit()

    def _fill_table(self):
        self.table.setRowCount(len(self.entries))
        for r, e in enumerate(self.entries):
            ref = self.ref_entry
            if ref is not None and ref.result is not None and e.result is not None:
                rm = f"{dist.rmsd(e.result.coords_std, ref.result.coords_std):.4f}"
            else:
                rm = "—"
            D = e.dist
            vals = [str(r + 1), e.dataset, "" if e.temperature is None else f"{e.temperature:g}", e.metal,
                    sym.pretty(e.symres.label), "" if e.result is None else f"{e.result.norm:.4f}",
                    f"{D.zeta:.4f}", dist.sci(D.delta), f"{D.sigma:.1f}", f"{D.theta:.1f}",
                    f"{D.theta_min:.1f}", f"{D.theta_max:.1f}", f"{D.volume:.3f}", rm, e.role]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if c == 0:
                    it.setForeground(SERIES_COLORS[r % len(SERIES_COLORS)])
                    f = it.font()
                    f.setBold(True)
                    it.setFont(f)
                if 5 <= c <= 13:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if c == 14 and e.role == "reference":
                    it.setForeground(QColor(0, 120, 0))
                self.table.setItem(r, c, it)

    def export_txt(self):
        shown = [e for e in self.entries if e.result is not None]
        if not shown:
            return
        p = export_q_dialog(self, [e.name for e in shown], [e.result.Q for e in shown], "Q_compare.txt")
        if p:
            self.msg.setText(self.msg.text() + f"<br>Saved: {p}")

    def copy_table(self):
        shown = [e for e in self.entries if e.result is not None]
        if not shown:
            return
        head = "Mode\tIrrep\t" + "\t".join(e.name for e in shown)
        rows = [head, "Symmetry\t\t" + "\t".join(e.symres.label for e in shown),
                "Frame\t\t" + "\t".join(e.role for e in shown)]
        for k in range(15):
            rows.append(f"Q{k + 1}\t{vv.IRREPS[k]}\t" + "\t".join(f"{e.result.Q[k]:.6f}" for e in shown))
        rows.append("|Q|\t(Q2–Q15)\t" + "\t".join(f"{e.result.norm:.6f}" for e in shown))
        QGuiApplication.clipboard().setText("\n".join(rows))
