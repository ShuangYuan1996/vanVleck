"""3D molecule view (pure QPainter rendering, no OpenGL required).

Controls: left-drag rotate · wheel / right-drag zoom · middle-drag or Shift+left-drag pan
          left-click an atom -> atomClicked(index); click empty space -> atomClicked(-1)
"""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import (QBrush, QColor, QFont, QLinearGradient, QPainter, QPen,
                           QPolygonF, QRadialGradient)
from PySide6.QtWidgets import QWidget

from .. import elements as el
from ..structure import Structure

STYLES = {
    "ball_stick": "Ball and stick",
    "stick": "Capped sticks",
    "wireframe": "Wireframe",
    "spacefill": "Spacefill",
}

HIGHLIGHT = QColor(255, 200, 0)


def _rot(axis: str, ang: float) -> np.ndarray:
    c, s = np.cos(ang), np.sin(ang)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


class MoleculeView(QWidget):
    atomClicked = Signal(int)
    atomHovered = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setMinimumSize(400, 300)
        self.setFocusPolicy(Qt.StrongFocus)
        self.structure: Structure | None = None
        self.style = "ball_stick"
        self.show_h = True
        self.show_labels = False
        self.show_cell = False
        self.show_polyhedron = True
        self.show_symmetry = True
        self.show_frame = True
        self.frame = None
        self.sym_elements = []
        self.dark_bg = False
        self.R = np.eye(3)
        self.center = np.zeros(3)
        self.scale = 30.0
        self.pan = np.zeros(2)
        self.sphere_idx: set[int] = set()
        self.metal_idx = -1
        self.faces: list[tuple[int, int, int]] = []
        self.hover_idx = -1
        self._press_pos = None
        self._press_btn = None
        self._dragged = False
        self._colors: dict[str, QColor] = {}

    # ------------------------------------------------------------------ API
    def set_structure(self, st: Structure | None):
        self.structure = st
        self._deg = [0] * (len(st.atoms) if st else 0)
        if st:
            for i, j in st.bonds:
                self._deg[i] += 1
                self._deg[j] += 1
        self.clear_highlight()
        self.reset_view()

    def set_style(self, style: str):
        self.style = style
        self.update()

    def set_option(self, name: str, value: bool):
        setattr(self, name, value)
        self.update()

    def set_highlight(self, metal_idx: int, ligand_idx: list[int], faces=None):
        self.metal_idx = metal_idx
        self.sphere_idx = {metal_idx, *ligand_idx}
        self.faces = list(faces or [])
        self.update()

    def set_symmetry_elements(self, elements):
        """elements: [(kind, order, axis, center)] with kind in C / S / σ."""
        self.sym_elements = list(elements)
        self.update()

    def set_frame(self, R, center):
        """Standard frame: the rows of R are the x, y, z axes expressed in structure coordinates."""
        self.frame = None if R is None else (np.array(R), np.array(center))
        self.update()

    def clear_highlight(self):
        self.frame = None
        self.sym_elements = []
        self.metal_idx = -1
        self.sphere_idx = set()
        self.faces = []
        self.update()

    def reset_view(self):
        self.R = np.eye(3)
        self.pan = np.zeros(2)
        st = self.structure
        if st is not None and st.atoms:
            xyz = st.coords
            self.center = (xyz.min(0) + xyz.max(0)) / 2
            ext = np.linalg.norm(xyz - self.center, axis=1).max() + 2.0
            self.scale = 0.48 * min(self.width(), self.height()) / ext
        self.update()

    def view_along(self, axis: str):
        """View along crystal axis a/b/c or Cartesian x/y/z."""
        st = self.structure
        if st is None:
            return
        if st.cell is not None and axis in "abc":
            M = st.cell.matrix
            d = M[:, "abc".index(axis)]
        else:
            d = np.eye(3)["xyz".index(axis) if axis in "xyz" else 2]
        z = d / np.linalg.norm(d)
        up = np.array([0, 1.0, 0]) if abs(z[1]) < 0.9 else np.array([1.0, 0, 0])
        x = np.cross(up, z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        self.R = np.array([x, y, z])
        self.update()

    # ------------------------------------------------------------ projection / picking
    def _color(self, sym: str) -> QColor:
        c = self._colors.get(sym)
        if c is None:
            c = self._colors[sym] = QColor(el.color(sym))
        return c

    def _visible(self) -> np.ndarray:
        st = self.structure
        vis = np.ones(len(st.atoms), bool)
        if not self.show_h:
            vis &= np.array([a.element != "H" for a in st.atoms])
        return vis

    def _project(self) -> np.ndarray:
        v = (self.structure.coords - self.center) @ self.R.T
        out = np.empty_like(v)
        out[:, 0] = self.width() / 2 + self.pan[0] + v[:, 0] * self.scale
        out[:, 1] = self.height() / 2 + self.pan[1] - v[:, 1] * self.scale
        out[:, 2] = v[:, 2]
        return out

    def _atom_radius(self, sym: str) -> float:
        if self.style == "spacefill":
            return el.vdw_radius(sym)
        if self.style == "stick":
            return 0.15
        if self.style == "wireframe":
            return 0.0
        if sym == "H":
            return 0.16
        if sym in el.METALS:
            return 0.42
        return 0.28

    def pick(self, x: float, y: float) -> int:
        st = self.structure
        if st is None or not st.atoms:
            return -1
        p = self._project()
        vis = self._visible()
        best, best_z = -1, -1e9
        for i, a in enumerate(st.atoms):
            if not vis[i]:
                continue
            r = max(self._atom_radius(a.element) * self.scale, 6.0)
            if (p[i, 0] - x) ** 2 + (p[i, 1] - y) ** 2 <= r * r and p[i, 2] > best_z:
                best, best_z = i, p[i, 2]
        return best

    # ---------------------------------------------------------------- drawing
    def paintEvent(self, _):
        qp = QPainter(self)
        qp.setRenderHint(QPainter.Antialiasing)
        bg = QColor(20, 22, 28) if self.dark_bg else QColor(255, 255, 255)
        fg = QColor(230, 230, 230) if self.dark_bg else QColor(40, 40, 40)
        qp.fillRect(self.rect(), bg)
        st = self.structure
        if st is None or not st.atoms:
            qp.setPen(QColor(140, 140, 140))
            qp.setFont(QFont("Segoe UI", 12))
            qp.drawText(self.rect(), Qt.AlignCenter,
                        "Open or drop a .cif / .xyz file\n\nClick a complex to highlight its inner coordination sphere")
            return

        P = self._project()
        vis = self._visible()
        highlight = bool(self.sphere_idx)
        style = self.style
        s = self.scale

        if self.show_cell and st.cell is not None:
            self._draw_cell(qp, fg)

        items = []
        for i in range(len(st.atoms)):
            if vis[i]:
                items.append((P[i, 2], 1, i, -1))
        if style != "spacefill":
            for i, j in st.bonds:
                if vis[i] and vis[j]:
                    zm = (P[i, 2] + P[j, 2]) / 2
                    items.append(((P[i, 2] + zm) / 2, 0, i, j))
                    items.append(((P[j, 2] + zm) / 2, 0, j, i))
        if highlight and self.show_polyhedron:
            for k, f in enumerate(self.faces):
                items.append((float(np.mean(P[list(f), 2])), 2, k, -2))
        sym_geo = self._symmetry_geometry() if (highlight and self.show_symmetry) else []
        for k, g in enumerate(sym_geo):
            items.append((g[0], 3, k, -3))
        items.sort(key=lambda t: (t[0], t[1]))

        bond_r = 0.15 if style in ("stick", "ball_stick") else 0.0
        for z, kind, i, j in items:
            if kind == 3:
                self._draw_sym(qp, sym_geo[i])
            elif kind == 2:
                self._draw_face(qp, P, self.faces[i])
            elif kind == 0:
                faded = highlight and not (i in self.sphere_idx and j in self.sphere_idx)
                self._draw_half_bond(qp, P, i, j, bond_r * s, faded)
            else:
                faded = highlight and i not in self.sphere_idx
                self._draw_atom(qp, P, i, faded)

        # labels
        qp.setFont(QFont("Arial", 9, QFont.Bold))
        for i, a in enumerate(st.atoms):
            if not vis[i]:
                continue
            in_sphere = i in self.sphere_idx
            if not (in_sphere or (self.show_labels and a.element != "H" and not highlight)):
                continue
            r = self._atom_radius(a.element) * s
            qp.setPen(QColor(200, 30, 30) if in_sphere and not self.dark_bg else
                      (QColor(255, 220, 80) if in_sphere else fg))
            qp.drawText(QPointF(P[i, 0] + r * 0.7 + 2, P[i, 1] - r * 0.7 - 2), a.label)

        # standard frame axes
        if highlight and self.show_frame and self.frame is not None:
            R, c = self.frame
            qp.setFont(QFont("Arial", 11, QFont.Bold))
            col = QColor(235, 235, 235) if self.dark_bg else QColor(0, 0, 0)   # black (white on a dark background)
            pen = QPen(col, 3.5, Qt.SolidLine, Qt.RoundCap)
            for k, name in enumerate("xyz"):
                x, y, _ = self._to_screen(np.array([c, c + 4.6 * R[k], c - 4.6 * R[k]]))
                qp.setPen(pen)
                qp.drawLine(QPointF(x[2], y[2]), QPointF(x[1], y[1]))   # one solid line from the negative to the positive end
                d = np.array([x[1] - x[0], y[1] - y[0]])
                nd = np.linalg.norm(d)
                if nd > 1:                                              # arrowhead at the positive end
                    d /= nd
                    tip = np.array([x[1], y[1]])
                    for s_ in (1, -1):
                        w = tip - 11 * d + s_ * 5 * np.array([-d[1], d[0]])
                        qp.drawLine(QPointF(*tip), QPointF(*w))
                qp.drawText(QPointF(x[1] + 6, y[1] + 4), name)

        # hover ring
        if 0 <= self.hover_idx < len(st.atoms) and vis[self.hover_idx]:
            i = self.hover_idx
            r = max(self._atom_radius(st.atoms[i].element) * s, 5) + 3
            qp.setBrush(Qt.NoBrush)
            qp.setPen(QPen(QColor(0, 160, 255), 2))
            qp.drawEllipse(QPointF(P[i, 0], P[i, 1]), r, r)

        if self.show_cell and st.cell is not None:
            self._draw_axes_triad(qp, fg)

    def _draw_atom(self, qp, P, i, faded):
        a = self.structure.atoms[i]
        x, y = P[i, 0], P[i, 1]
        col = self._color(a.element)
        alpha = 45 if faded else 255
        if self.style == "wireframe":
            r = 3.0 if not self._deg[i] else 0.0
            if i in self.sphere_idx:
                r = 5.0
            if r:
                c = QColor(col)
                c.setAlpha(alpha)
                qp.setPen(Qt.NoPen)
                qp.setBrush(c)
                qp.drawEllipse(QPointF(x, y), r, r)
            if i in self.sphere_idx:
                qp.setBrush(Qt.NoBrush)
                qp.setPen(QPen(HIGHLIGHT, 2))
                qp.drawEllipse(QPointF(x, y), r + 3, r + 3)
            return
        r = self._atom_radius(a.element) * self.scale
        if i in self.sphere_idx:
            r *= 1.15
        g = QRadialGradient(QPointF(x - r * 0.35, y - r * 0.35), r * 1.3)
        c0, c1, c2 = QColor(col).lighter(190), QColor(col), QColor(col).darker(230)
        for c in (c0, c1, c2):
            c.setAlpha(alpha)
        g.setColorAt(0.0, c0)
        g.setColorAt(0.45, c1)
        g.setColorAt(1.0, c2)
        qp.setBrush(QBrush(g))
        edge = QColor(60, 60, 60, 120 if not faded else 25)
        qp.setPen(QPen(edge, 0.6))
        qp.drawEllipse(QPointF(x, y), r, r)
        if i in self.sphere_idx:
            qp.setBrush(Qt.NoBrush)
            qp.setPen(QPen(HIGHLIGHT, 2.5))
            qp.drawEllipse(QPointF(x, y), r + 2.5, r + 2.5)

    def _draw_half_bond(self, qp, P, i, j, width, faded):
        a = self.structure.atoms[i]
        x1, y1 = P[i, 0], P[i, 1]
        xm, ym = (P[i, 0] + P[j, 0]) / 2, (P[i, 1] + P[j, 1]) / 2
        col = QColor(self._color(a.element))
        is_ml = self.metal_idx >= 0 and self.metal_idx in (i, j) and i in self.sphere_idx and j in self.sphere_idx
        alpha = 40 if faded else 255
        if self.style == "wireframe" or width < 1.0:
            col.setAlpha(alpha)
            qp.setPen(QPen(col, 3.0 if is_ml else 1.6, Qt.SolidLine, Qt.RoundCap))
            qp.drawLine(QPointF(x1, y1), QPointF(xm, ym))
            return
        if is_ml:
            width *= 1.3
        dx, dy = xm - x1, ym - y1
        L = (dx * dx + dy * dy) ** 0.5
        if L < 0.5:
            return
        nx, ny = -dy / L * width / 2, dx / L * width / 2
        g = QLinearGradient(QPointF(xm + nx, ym + ny), QPointF(xm - nx, ym - ny))
        cs = [col.darker(200), col.lighter(150), col.darker(200)]
        for c in cs:
            c.setAlpha(alpha)
        g.setColorAt(0.0, cs[0])
        g.setColorAt(0.4, cs[1])
        g.setColorAt(1.0, cs[2])
        cap = Qt.RoundCap if self.style == "stick" else Qt.FlatCap
        qp.setPen(QPen(QBrush(g), width, Qt.SolidLine, cap))
        qp.drawLine(QPointF(x1, y1), QPointF(xm, ym))

    # ---- symmetry elements ----
    _AXIS_COLORS = {2: QColor(30, 120, 255), 3: QColor(0, 170, 90), 4: QColor(220, 30, 160), 6: QColor(130, 60, 200)}

    def _to_screen(self, pts):
        v = (np.atleast_2d(pts) - self.center) @ self.R.T
        x = self.width() / 2 + self.pan[0] + v[:, 0] * self.scale
        y = self.height() / 2 + self.pan[1] - v[:, 1] * self.scale
        return x, y, v[:, 2]

    def _symmetry_geometry(self):
        """Turn symmetry elements into drawing primitives (depth, type, data); one label per axis."""
        geo, axes = [], []          # axes: [(axis vector, {symbols})]
        L = 3.4
        for kind, order, axis, c in self.sym_elements:
            if kind in ("C", "S"):
                for a_, labs in axes:
                    if abs(abs(np.dot(axis, a_)) - 1) < 1e-3:
                        labs.add((kind, order))
                        break
                else:
                    axes.append((axis, {(kind, order)}))
            elif kind == "σ":
                n = axis
                u = np.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0])
                u /= np.linalg.norm(u)
                w = np.cross(n, u)
                ring = [c + 2.9 * (np.cos(t) * u + np.sin(t) * w) for t in np.linspace(0, 2 * np.pi, 40, endpoint=False)]
                x, y, z = self._to_screen(np.array(ring))
                geo.append((float(self._to_screen(c)[2][0]) - 0.02, "plane", (x, y)))
        axes = [(a_, "/".join(f"{k}{o}" for k, o in sorted(labs, key=lambda t: (t[0] != "C", -t[1]))
                              if not (k == "C" and any(k2 == "C" and o2 > o and o2 % o == 0 for k2, o2 in labs))))
                for a_, labs in axes]
        for a, lab in axes:
            c = self.sym_elements[0][3]
            x, y, z = self._to_screen(np.array([c - L * a, c + L * a]))
            order = max(int(t[1:]) for t in lab.split("/") if t[0] == "C") if "C" in lab else 4
            geo.append((float(z.max()) + 0.01, "axis", (x, y, z, lab, order)))
        return geo

    def _draw_sym(self, qp, g):
        _, typ, data = g
        if typ == "plane":
            x, y = data
            c = QColor(120, 120, 255, 45)
            qp.setBrush(c)
            qp.setPen(QPen(QColor(90, 90, 220, 130), 1.0))
            qp.drawPolygon(QPolygonF([QPointF(a, b) for a, b in zip(x, y)]))
        else:
            x, y, z, lab, order = data
            col = self._AXIS_COLORS.get(order, QColor(80, 80, 80))
            qp.setPen(QPen(col, 2.2, Qt.DashLine, Qt.RoundCap))
            qp.drawLine(QPointF(x[0], y[0]), QPointF(x[1], y[1]))
            k = int(np.argmax(z))
            qp.setFont(QFont("Arial", 9, QFont.Bold))
            qp.drawText(QPointF(x[k] + 4, y[k] - 4), lab)

    def _draw_face(self, qp, P, f):
        poly = QPolygonF([QPointF(P[k, 0], P[k, 1]) for k in f])
        c = QColor(self._color(self.structure.atoms[self.metal_idx].element))
        c.setAlpha(55)
        qp.setBrush(c)
        e = QColor(c).darker(150)
        e.setAlpha(140)
        qp.setPen(QPen(e, 1.0))
        qp.drawPolygon(poly)

    def _cell_corners(self):
        M = self.structure.cell.matrix
        fr = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float)
        v = (fr @ M.T - self.center) @ self.R.T
        x = self.width() / 2 + self.pan[0] + v[:, 0] * self.scale
        y = self.height() / 2 + self.pan[1] - v[:, 1] * self.scale
        return fr, x, y

    def _draw_cell(self, qp, fg):
        fr, x, y = self._cell_corners()
        axis_col = {0: QColor(220, 40, 40), 1: QColor(40, 170, 40), 2: QColor(40, 80, 230)}
        for a in range(8):
            for b in range(a + 1, 8):
                diff = fr[b] - fr[a]
                if np.abs(diff).sum() != 1:
                    continue
                ax = int(np.argmax(diff))
                if a == 0:
                    pen = QPen(axis_col[ax], 2)
                else:
                    c = QColor(fg)
                    c.setAlpha(110)
                    pen = QPen(c, 1)
                qp.setPen(pen)
                qp.drawLine(QPointF(x[a], y[a]), QPointF(x[b], y[b]))

    def _draw_axes_triad(self, qp, fg):
        M = self.structure.cell.matrix
        ox, oy, L = 45, self.height() - 45, 30
        qp.setFont(QFont("Arial", 10, QFont.Bold))
        cols = [QColor(220, 40, 40), QColor(40, 170, 40), QColor(40, 80, 230)]
        for k, name in enumerate("abc"):
            d = self.R @ (M[:, k] / np.linalg.norm(M[:, k]))
            ex, ey = ox + d[0] * L, oy - d[1] * L
            qp.setPen(QPen(cols[k], 2.5))
            qp.drawLine(QPointF(ox, oy), QPointF(ex, ey))
            qp.drawText(QPointF(ex + d[0] * 8 - 4, ey - d[1] * 8 + 4), name)

    # ---------------------------------------------------------------- mouse
    def mousePressEvent(self, e):
        self._press_pos = e.position()
        self._last = e.position()
        self._press_btn = e.button()
        self._dragged = False

    def mouseMoveEvent(self, e):
        pos = e.position()
        if self._press_btn is None:
            idx = self.pick(pos.x(), pos.y()) if self.structure else -1
            if idx != self.hover_idx:
                self.hover_idx = idx
                self.atomHovered.emit(idx)
                self.update()
            return
        dx, dy = pos.x() - self._last.x(), pos.y() - self._last.y()
        self._last = pos
        if (pos - self._press_pos).manhattanLength() > 3:
            self._dragged = True
        btn = self._press_btn
        if btn == Qt.LeftButton and not (e.modifiers() & Qt.ShiftModifier):
            self.R = _rot("y", dx * 0.01) @ _rot("x", dy * 0.01) @ self.R
        elif btn == Qt.MiddleButton or btn == Qt.LeftButton:
            self.pan += (dx, dy)
        elif btn == Qt.RightButton:
            self.scale *= float(np.exp(-dy * 0.01))
        self.update()

    def mouseReleaseEvent(self, e):
        if e.button() == Qt.LeftButton and not self._dragged and self.structure:
            self.atomClicked.emit(self.pick(e.position().x(), e.position().y()))
        self._press_btn = None

    def wheelEvent(self, e):
        f = 1.0 + 0.1 * np.sign(e.angleDelta().y())
        self.scale *= f
        self.update()

    def leaveEvent(self, _):
        if self.hover_idx != -1:
            self.hover_idx = -1
            self.atomHovered.emit(-1)
            self.update()
