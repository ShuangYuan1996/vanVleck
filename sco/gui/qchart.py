"""Bar charts of the van Vleck coordinates Q1–Q15 (drawn with QPainter, no extra dependencies)."""
from __future__ import annotations

import re

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QGuiApplication, QPainter, QPen
from PySide6.QtWidgets import (QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout,
                               QWidget)

# colour by irreducible representation
IRREP_COLORS = {
    "A1g": QColor(120, 120, 120), "Eg": QColor(142, 68, 173), "T2g": QColor(230, 126, 34),
    "T2u": QColor(41, 98, 214), "T1u(I)": QColor(39, 160, 80), "T1u(II)": QColor(22, 120, 110),
}


def fmt3(v: float, zero_tol: float) -> str:
    """Three significant figures; "0" when |v| is below the tolerance."""
    if abs(v) < zero_tol:
        return "0"
    s = f"{v:.3g}"
    if "e" in s:                      # scientific notation for very small / large numbers
        m, e = s.split("e")
        s = f"{m}×10^{int(e)}"
    return s


SERIES_COLORS = [QColor(c) for c in ("#2a62d6", "#e0662b", "#2ca05a", "#8e44ad", "#c0392b", "#16a0a0",
                                        "#b8860b", "#555555", "#e84393", "#1e8449")]


class QBarChart(QWidget):
    """Q1–Q15 bar chart: one series is coloured by irrep; several series are grouped, coloured per dataset, with a legend."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(760, 420)
        self.Qs = None
        self.names = []
        self.irreps = []
        self.zero_tol = 1e-5
        self.title = ""
        self.subtitle = ""

    def set_data(self, Q, irreps, zero_tol, title="", subtitle=""):
        self.set_series([Q], [""], irreps, zero_tol, title, subtitle)

    def set_series(self, Qs, names, irreps, zero_tol, title="", subtitle=""):
        self.subtitle = subtitle
        self.Qs = [np.asarray(q, float) for q in Qs] if len(Qs) else None
        self.names = list(names)
        self.irreps = list(irreps)
        self.zero_tol = zero_tol
        self.title = title
        self.update()

    def paintEvent(self, _):
        qp = QPainter(self)
        qp.setRenderHint(QPainter.Antialiasing)
        qp.fillRect(self.rect(), QColor(255, 255, 255))
        if not self.Qs:
            qp.setPen(QColor(150, 150, 150))
            qp.drawText(self.rect(), Qt.AlignCenter, "No data")
            return
        multi = len(self.Qs) > 1
        W, H = self.width(), self.height()
        two_line = any("\n" in n for n in self.names)
        ncol = 2 if two_line else 3
        row_h = 36 if two_line else 22
        legend_h = row_h * ((len(self.names) + ncol - 1) // ncol) if multi else 0
        sub_lines = self.subtitle.split("\n") if self.subtitle else []
        sub_h = 20 * len(sub_lines)
        left, right, top, bottom = 70, 20, 50 + legend_h + sub_h, 60
        pw, ph = W - left - right, H - top - bottom

        allQ = np.concatenate(self.Qs)
        vmax = max(float(allQ.max()), 0.0)
        vmin = min(float(allQ.min()), 0.0)
        span = vmax - vmin if vmax > vmin else 1.0
        pad = 0.18 if multi else 0.12           # room for the value labels
        vmax += pad * span
        vmin -= pad * span if vmin < 0 else 0.0
        y_of = lambda v: top + (vmax - v) / (vmax - vmin) * ph

        qp.setPen(QColor(30, 30, 30))
        qp.setFont(QFont("Segoe UI", 11, QFont.Bold))
        qp.drawText(QRectF(0, 8, W, 26), Qt.AlignHCenter, self.title)
        if self.subtitle:
            qp.setPen(QColor(70, 70, 70))
            qp.setFont(QFont("Segoe UI", 10))
            for li, line in enumerate(sub_lines):
                qp.drawText(QRectF(0, 32 + 20 * li, W, 20), Qt.AlignHCenter, line)

        # legend
        if multi:
            qp.setFont(QFont("Segoe UI", 9))
            colw = (W - 2 * left) / ncol
            for i, name in enumerate(self.names):
                r_, c_ = divmod(i, ncol)
                x0, y0_ = left + c_ * colw, 38 + sub_h + r_ * row_h
                qp.setPen(Qt.NoPen)
                qp.setBrush(SERIES_COLORS[i % len(SERIES_COLORS)])
                qp.drawRect(QRectF(x0, y0_ + 3, 14, 12))
                first, _, second = name.partition("\n")
                qp.setPen(QColor(30, 30, 30))
                qp.setFont(QFont("Segoe UI", 9))
                qp.drawText(QRectF(x0 + 20, y0_, colw - 24, 18), Qt.AlignLeft | Qt.AlignVCenter, first)
                if second:
                    qp.setPen(QColor(90, 90, 90))
                    qp.setFont(QFont("Segoe UI", 8))
                    qp.drawText(QRectF(x0 + 20, y0_ + 16, colw - 24, 16), Qt.AlignLeft | Qt.AlignVCenter, second)

        # axes and ticks
        qp.setFont(QFont("Segoe UI", 8))
        fm = QFontMetricsF(qp.font())
        for t in _nice_ticks(vmin, vmax):
            y = y_of(t)
            qp.setPen(QPen(QColor(232, 232, 232), 1))
            qp.drawLine(QPointF(left, y), QPointF(left + pw, y))
            qp.setPen(QColor(80, 80, 80))
            lab = f"{t:g}"
            qp.drawText(QPointF(left - 6 - fm.horizontalAdvance(lab), y + 4), lab)
        qp.setPen(QPen(QColor(60, 60, 60), 1.2))
        qp.drawLine(QPointF(left, top), QPointF(left, top + ph))
        y0 = y_of(0.0)
        qp.drawLine(QPointF(left, y0), QPointF(left + pw, y0))
        qp.save()
        qp.translate(18, top + ph / 2)
        qp.rotate(-90)
        qp.setFont(QFont("Segoe UI", 9))
        qp.drawText(QRectF(-60, -10, 120, 20), Qt.AlignCenter, "Q (Å)")
        qp.restore()

        n = len(self.Qs[0])
        ns = len(self.Qs)
        slot = pw / n
        group_w = slot * (0.62 if ns == 1 else 0.84)
        bw = group_w / ns
        vertical = ns >= 3 or bw < 30               # vertical value labels when bars are narrow
        for k in range(n):
            cx = left + slot * (k + 0.5)
            irr = self.irreps[k] if k < len(self.irreps) else ""
            for si, Q in enumerate(self.Qs):
                v = Q[k]
                col = SERIES_COLORS[si % len(SERIES_COLORS)] if multi else IRREP_COLORS.get(irr, QColor(100, 100, 100))
                bx = cx - group_w / 2 + si * bw
                y1 = y_of(v)
                nonzero = abs(v) >= self.zero_tol
                if nonzero:
                    qp.setPen(Qt.NoPen)
                    qp.setBrush(col)
                    qp.drawRect(QRectF(bx + 0.5, min(y0, y1), bw - 1, max(abs(y1 - y0), 1.0)))
                lab = fmt3(v, self.zero_tol)
                qp.setFont(QFont("Segoe UI", 7 if multi else 8, QFont.Bold))
                fm2 = QFontMetricsF(qp.font())
                tw = fm2.horizontalAdvance(lab)
                qp.setPen(QColor(30, 30, 30) if nonzero else QColor(150, 150, 150))
                up = v >= 0 or not nonzero
                yb = y0 if not nonzero else y1
                if vertical:
                    qp.save()
                    qp.translate(bx + bw / 2 + 3, yb - 3 if up else yb + 3)
                    qp.rotate(-90)
                    qp.drawText(QPointF(0 if up else -tw, 0), lab)
                    qp.restore()
                else:
                    qp.drawText(QPointF(bx + bw / 2 - tw / 2, (yb - 5) if up else (yb + 13)), lab)
            # with exactly two series, a Q that is zero in one and non-zero in the other is marked red bold
            zeros = [abs(Q[k]) < self.zero_tol for Q in self.Qs]
            mixed = len(self.Qs) == 2 and zeros[0] != zeros[1]
            qp.setFont(QFont("Segoe UI", 10 if mixed else 9, QFont.Bold if mixed else QFont.Normal))
            qp.setPen(QColor(210, 0, 0) if mixed else QColor(40, 40, 40))
            qp.drawText(QRectF(cx - slot / 2, top + ph + 6, slot, 16), Qt.AlignHCenter, f"Q{k + 1}")
            qp.setFont(QFont("Segoe UI", 7))
            qp.setPen(IRREP_COLORS.get(irr, QColor(100, 100, 100)))
            qp.drawText(QRectF(cx - slot / 2, top + ph + 22, slot, 14), Qt.AlignHCenter, irr)


def _nice_ticks(lo, hi, n=6):
    span = hi - lo
    if span <= 0:
        return [0.0]
    raw = span / n
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=raw)
    start = np.ceil(lo / step) * step
    out, t = [], start
    while t <= hi + 1e-12:
        out.append(0.0 if abs(t) < step * 1e-9 else float(t))
        t += step
    return out


def write_q_txt(path: str, names, Qs) -> None:
    """Tab-separated text: header Name, Q1…Q15, |Q|; then one row per molecule (name first).
    |Q| = sqrt(Q2² + … + Q15²)."""
    lines = ["Name\t" + "\t".join(f"Q{k + 1}" for k in range(15)) + "\t|Q|"]
    for name, Q in zip(names, Qs):
        clean = str(name).replace("\t", " ")
        nrm = float(np.linalg.norm(np.asarray(Q, float)[1:15]))
        lines.append(clean + "\t" + "\t".join(f"{(0.0 if abs(v) < 5e-7 else float(v)):.6f}" for v in Q)
                     + f"\t{nrm:.6f}")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def export_q_dialog(parent, names, Qs, default_name="Q.txt") -> str | None:
    """Ask for a file name and write the .txt; return the path (None if cancelled)."""
    if not len(Qs):
        return None
    default_name = re.sub(r'[\\/:*?"<>|]+', "_", default_name)
    path, _ = QFileDialog.getSaveFileName(parent, "Export Q values", default_name, "Text files (*.txt)")
    if not path:
        return None
    if not path.lower().endswith(".txt"):
        path += ".txt"
    try:
        write_q_txt(path, names, Qs)
    except OSError as exc:
        QMessageBox.critical(parent, "Export failed", f"Could not write file:\n{path}\n\n{exc}")
        return None
    return path


class QChartDialog(QDialog):
    """Modeless pop-up, refreshed each time a coordination sphere is selected."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("van Vleck coordinates Q1–Q15")
        self.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        self.resize(900, 520)
        lay = QVBoxLayout(self)
        self.chart = QBarChart()
        lay.addWidget(self.chart, 1)
        row = QHBoxLayout()
        self.info = QLabel("")
        row.addWidget(self.info, 1)
        btn_exp = QPushButton("Export .txt")
        btn_exp.setToolTip("Save Q1–Q15 as a tab-separated .txt file")
        btn_exp.clicked.connect(self.export_txt)
        row.addWidget(btn_exp)
        btn_copy = QPushButton("Copy image")
        btn_copy.clicked.connect(lambda: QGuiApplication.clipboard().setPixmap(self.chart.grab()))
        row.addWidget(btn_copy)
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        row.addWidget(btn_close)
        lay.addLayout(row)

    def show_result(self, Q, irreps, zero_tol, title, info="", name="", subtitle=""):
        self.name = name or title
        self.chart.set_data(Q, irreps, zero_tol, title, subtitle)
        self.info.setText(info)
        self.show()
        self.raise_()

    def export_txt(self):
        if self.chart.Qs:
            p = export_q_dialog(self, [self.name], self.chart.Qs, f"Q_{self.name}.txt")
            if p:
                self.info.setText(f"Saved: {p}")
