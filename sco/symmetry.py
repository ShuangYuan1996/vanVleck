"""Point-group analysis of an ML6 inner coordination sphere.

Every symmetry operation must leave the unique metal atom fixed, so the metal is the origin.
Candidate axes are generated from the donor directions (single vectors, pairwise sums,
differences and cross products, sums of triples). C2/C3/C4/C6 rotations, mirrors, the inversion
and S4/S6 are tested on each axis; an operation is accepted when it maps the 7 atoms onto atoms
of the same element within `tol`. The accepted operations are closed into a group and named.

D2h, C2v, C2h, C2 and Cs are split into variants according to the number of donor atoms lying
on their symmetry elements: D2h^v / D2h^d, C2v^v / C2v^d (/ C2v^m), C2^z / C2^xy,
Cs^z / Cs^xy (written Cs^{x-ybar}), and C2h / C2h^xy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np

DEFAULT_TOL = 1e-5   # Å


# ----------------------------------------------------------------------------
# basic geometry
# ----------------------------------------------------------------------------
def rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    n = axis / np.linalg.norm(axis)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def reflection(normal: np.ndarray) -> np.ndarray:
    n = normal / np.linalg.norm(normal)
    return np.eye(3) - 2.0 * np.outer(n, n)


def classify_op(R: np.ndarray) -> tuple[str, int, np.ndarray | None]:
    """Return (kind, order, axis) with kind in E, C, i, σ, S."""
    det = np.linalg.det(R)
    if det > 0:
        c = np.clip((np.trace(R) - 1) / 2, -1, 1)
        ang = np.arccos(c)
        if ang < 1e-3:
            return "E", 1, None
        w, v = np.linalg.eig(R)
        axis = np.real(v[:, np.argmin(np.abs(w - 1))])
        return "C", int(round(2 * np.pi / ang)), axis / np.linalg.norm(axis)
    P = -R                                      # improper rotation = inversion × proper rotation
    c = np.clip((np.trace(P) - 1) / 2, -1, 1)
    ang = np.arccos(c)                          # rotation angle of P
    if ang < 1e-3:
        return "i", 2, None
    w, v = np.linalg.eig(R)
    axis = np.real(v[:, np.argmin(np.abs(w + 1))])
    axis /= np.linalg.norm(axis)
    if abs(ang - np.pi) < 1e-3:
        return "σ", 2, axis                     # mirror; axis is its normal
    # S_n: rotation angle theta = pi - ang
    theta = np.pi - ang
    return "S", int(round(2 * np.pi / theta)), axis


# ----------------------------------------------------------------------------
# result data structures
# ----------------------------------------------------------------------------
@dataclass
class SymOp:
    matrix: np.ndarray
    kind: str            # E, C, i, σ, S
    order: int
    axis: np.ndarray | None
    perm: list[int]      # atom i -> perm[i]
    deviation: float     # largest atomic deviation (Å)

    @property
    def fixed_ligands(self) -> list[int]:
        """Donor atoms (indices 1–6) lying on this symmetry element."""
        return [i for i in range(1, len(self.perm)) if self.perm[i] == i]

    @property
    def symbol(self) -> str:
        if self.kind == "E":
            return "E"
        if self.kind == "i":
            return "i"
        if self.kind == "σ":
            return "σ"
        return f"{self.kind}{self.order}"


@dataclass
class SymmetryResult:
    point_group: str                 # Schoenflies symbol, e.g. 'D2d'
    label: str                       # symbol with variant superscript, e.g. 'C2v^d'
    order: int
    ops: list[SymOp]
    tol: float
    max_deviation: float             # largest deviation among the accepted operations
    standard: bool                   # one of the 19 standard ML6 subgroups?
    notes: list[str] = field(default_factory=list)

    def elements_summary(self) -> list[str]:
        """List of symmetry elements (powers on one axis merged), for display."""
        out, seen = [], []
        for op in sorted(self.ops, key=lambda o: ({"E": 0, "C": 1, "S": 2, "σ": 3, "i": 4}[o.kind], -o.order)):
            if op.kind in ("E", "i"):
                out.append(f"{op.symbol}   (dev {op.deviation:.4f} Å)")
                continue
            key = (op.kind, op.order)
            if any(k == key and abs(abs(np.dot(a, op.axis)) - 1) < 1e-3 for k, a in seen):
                continue
            # skip lower powers on an axis already listed (e.g. C2 on a C4 axis)
            seen.append((key, op.axis))
            n = len(op.fixed_ligands)
            where = "normal" if op.kind == "σ" else "axis"
            out.append(f"{op.symbol}   {where} ({op.axis[0]:+.3f}, {op.axis[1]:+.3f}, {op.axis[2]:+.3f}),"
                       f" {n} ligand(s) on element, dev {op.deviation:.4f} Å")
        return out


# the 19 standard ML6 subgroups (with variants)
STANDARD_GROUPS = [
    "Oh", "D4h", "D3d", "D2d", "C4v", "D2h^d", "D2h^v", "C3v", "D2", "C2v^v",
    "C2v^d", "C2h", "C3", "C2^xy", "C2^z", "Cs^xy", "Cs^z", "Ci", "C1",
]


# ----------------------------------------------------------------------------
# core algorithm
# ----------------------------------------------------------------------------
def _match(Y: np.ndarray, X: np.ndarray, elems: list[str]) -> tuple[list[int] | None, float]:
    """Match transformed coordinates Y to X (same element, nearest neighbour); return permutation and max deviation."""
    n = len(X)
    perm, used, dev = [-1] * n, set(), 0.0
    for i in range(n):
        best, bd = -1, 1e9
        for j in range(n):
            if elems[j] != elems[i] or j in used:
                continue
            d = np.linalg.norm(Y[i] - X[j])
            if d < bd:
                best, bd = j, d
        if best < 0:
            return None, 1e9
        perm[i] = best
        used.add(best)
        dev = max(dev, bd)
    return perm, dev


def _candidate_axes(rel: np.ndarray) -> list[np.ndarray]:
    u = [v / np.linalg.norm(v) for v in rel[1:] if np.linalg.norm(v) > 1e-6]
    cands = list(u)
    for a, b in combinations(u, 2):
        cands += [a + b, a - b, np.cross(a, b)]
    for a, b, c in combinations(u, 3):
        cands.append(a + b + c)
    for (a, b), (c, d) in combinations(list(combinations(u, 2)), 2):
        cands.append(np.cross(a - b, c - d))
    out = []
    for v in cands:
        nv = np.linalg.norm(v)
        if nv < 1e-3:
            continue
        v = v / nv
        if not any(abs(np.dot(v, w)) > 1 - 1e-9 for w in out):
            out.append(v)
    return out


def _op_matrices(axis: np.ndarray):
    """All operation matrices to test about one axis."""
    sig = reflection(axis)
    for k in (2, 3, 4, 6):
        for m in range(1, k):
            if np.gcd(m, k) != 1:
                continue
            R = rotation(axis, 2 * np.pi * m / k)
            yield R
            if k in (4, 6, 3):
                yield sig @ R       # S4, S6, S3
    yield sig


def _same(A: np.ndarray, B: np.ndarray) -> bool:
    return np.linalg.norm(A - B) < 0.25


def find_operations(coords: np.ndarray, elems: list[str], tol: float = DEFAULT_TOL) -> list[SymOp]:
    X = np.asarray(coords, float)
    rel = X - X[0]                                # metal at the origin
    found: list[SymOp] = []

    def test(R):
        perm, dev = _match(rel @ R.T, rel, elems)
        if perm is None or dev > tol or perm[0] != 0:
            return
        for k, f in enumerate(found):
            if _same(R, f.matrix):
                if dev < f.deviation:      # keep the one with the smaller deviation (more accurate axis)
                    kind, order, axis = classify_op(R)
                    found[k] = SymOp(R, kind, order, axis, perm, dev)
                return
        kind, order, axis = classify_op(R)
        found.append(SymOp(R, kind, order, axis, perm, dev))

    test(np.eye(3))
    test(-np.eye(3))
    for ax in _candidate_axes(rel):
        for R in _op_matrices(ax):
            test(R)

    # group closure: test any product not yet in the set
    changed = True
    while changed and len(found) <= 48:
        changed = False
        for a in list(found):
            for b in list(found):
                P = a.matrix @ b.matrix
                if not any(_same(P, f.matrix) for f in found):
                    n0 = len(found)
                    test(P)
                    changed |= len(found) > n0

    # at large tolerance the set may not close: drop the worst operation until it is a group
    def closed(ops):
        return all(any(_same(a.matrix @ b.matrix, c.matrix) for c in ops) for a in ops for b in ops)

    while not closed(found):
        worst = max((o for o in found if o.kind != "E"), key=lambda o: o.deviation)
        found = [o for o in found if o is not worst]
    return found


def _schoenflies(ops: list[SymOp]) -> str:
    n = len(ops)
    cnt = lambda kind, order=None: sum(1 for o in ops if o.kind == kind and (order is None or o.order == order))
    inv = cnt("i") > 0
    nC2, nC3, nC4, nC6 = cnt("C", 2), cnt("C", 3), cnt("C", 4), cnt("C", 6)
    nσ, nS4, nS6, nS3 = cnt("σ"), cnt("S", 4), cnt("S", 6), cnt("S", 3)
    if n == 1:
        return "C1"
    if n == 2:
        return "Ci" if inv else ("Cs" if nσ else "C2")
    if n == 3:
        return "C3"
    if n == 4:
        if nC4:
            return "C4"
        if nS4:
            return "S4"
        if inv:
            return "C2h"
        return "C2v" if nσ == 2 else "D2"
    if n == 6:
        if inv:
            return "C3i"
        if nσ == 3:
            return "C3v"
        if nS3:
            return "C3h"
        if nC6:
            return "C6"
        return "D3"
    if n == 8:
        if inv:
            return "C4h" if nC4 else "D2h"
        if nC4:
            return "C4v" if nσ else "D4"
        return "D2d"
    if n == 12:
        if nC3 == 8:
            return "Th" if inv else "T"
        if inv:
            return "D3d" if nC6 == 0 else "C6h"
        if nσ == 6:
            return "C6v"
        if nS3:
            return "D3h"
        return "D6"
    if n == 16:
        return "D4h"
    if n == 24:
        if nC3 == 8:
            if inv:
                return "Th" if not nC4 else "?"
            return "O" if nC4 else "Td"
        if inv:
            return "D6h"
        return "?"
    if n == 48:
        return "Oh"
    return f"?(order {n})"


def _variant(pg: str, ops: list[SymOp]) -> tuple[str, list[str]]:
    """Add the variant superscript."""
    notes = []
    C2 = [o for o in ops if o.kind == "C" and o.order == 2]
    σ = [o for o in ops if o.kind == "σ"]
    if pg == "D2h":
        on = sum(1 for o in C2 if len(o.fixed_ligands) == 2)
        return ("D2h^v" if on == 3 else "D2h^d"), notes
    if pg == "C2v":
        n_in = sorted(len(o.fixed_ligands) for o in σ)
        if n_in == [4, 4]:
            return "C2v^v", notes
        if n_in == [2, 2]:
            return "C2v^d", notes
        notes.append("C2 axis contains no ligand atoms; one mirror contains 4 ligands and the other 2 "
                     "(e.g. cis-MA4B2 type). This C2v variant is outside the 19 standard groups.")
        return "C2v^m", notes
    if pg == "C2":
        return ("C2^z" if len(C2[0].fixed_ligands) == 2 else "C2^xy"), notes
    if pg == "Cs":
        return ("Cs^z" if len(σ[0].fixed_ligands) >= 4 else "Cs^xy"), notes
    if pg == "C2h":
        if len(C2[0].fixed_ligands) == 2:
            return "C2h", notes
        notes.append("C2 axis contains no ligand atoms (C2h^xy); this variant is outside the 19 standard groups.")
        return "C2h^xy", notes
    if pg == "D3":
        notes.append("D3 (e.g. tris-chelate twist) is outside the 19 standard groups.")
    return pg, notes


def analyze(coords: np.ndarray, elems: list[str], tol: float = DEFAULT_TOL) -> SymmetryResult:
    """coords: (7, 3) with the metal in row 0; elems: the 7 element symbols."""
    ops = find_operations(coords, elems, tol)
    pg = _schoenflies(ops)
    label, notes = _variant(pg, ops)
    if len(set(elems[1:])) > 1:
        notes.append("Donor atoms are not all the same element; only operations mapping like atoms onto like atoms are kept.")
    return SymmetryResult(point_group=pg, label=label, order=len(ops), ops=ops, tol=tol,
                          max_deviation=max(o.deviation for o in ops), standard=label in STANDARD_GROUPS,
                          notes=notes)


def symmetry_vs_tolerance(coords, elems, tols=(1e-5, 1e-3, 0.01, 0.05, 0.1, 0.2)) -> list[tuple[float, str]]:
    """Point group at several tolerances, to judge pseudo-symmetry."""
    return [(t, analyze(coords, elems, t).label) for t in tols]


def pretty(label: str) -> str:
    """Unicode display form, e.g. 'C2v^d' -> 'C₂ᵥᵈ'."""
    sub = str.maketrans("0123456789hvsi", "₀₁₂₃₄₅₆₇₈₉ₕᵥₛᵢ")
    base, _, sup = label.partition("^")
    head, rest = base[0], base[1:]
    rest = rest.replace("d", "\u2093").translate(sub).replace("\u2093", "d")
    sups = {"xy": "ˣʸ", "z": "ᶻ", "d": "ᵈ", "v": "ᵛ", "m": "ᵐ"}
    if base == "Cs" and sup == "xy":
        return head + rest + "ˣʸ\u0304"          # Cs^{xȳ}
    return head + rest + sups.get(sup, f"^{sup}" if sup else "")


def html(label: str) -> str:
    """HTML form, e.g. 'C2v^d' -> 'C<sub>2v</sub><sup>d</sup>'."""
    base, _, sup = label.partition("^")
    if base == "Cs" and sup == "xy":
        sup = "x<span style='text-decoration:overline'>y</span>"
    return f"{base[0]}<sub>{base[1:]}</sub>" + (f"<sup>{sup}</sup>" if sup else "")
