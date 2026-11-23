"""Van Vleck coordinates Q1–Q15.

1. Each point group has a representative subgroup H_std of Oh in a standard orientation
   (metal at the origin):
   - a C4 or S4 axis is z;
   - a C3 axis lies along [111];
   - a single C2 axis through donors is z; one through no donor lies along [110];
   - a single mirror containing 4 donors is the xy plane (normal z); one containing 2 donors is
     the x = y plane (normal [1-10]);
   - an inversion centre is on the metal.
2. A proper rotation R is found such that R g R^T is in H_std for every symmetry operation g of
   the molecule (symmetry elements exactly in standard positions) and, under that constraint,
   the six donor directions are best aligned with ±x, ±y, ±z (a constrained Procrustes problem
   with an analytic solution).
3. Donors are numbered 1:+x, 2:+y, 3:+z, 4:−x, 5:−y, 6:−z and inserted into the definitions.
4. Group theory (the invariant subspace of H_std) gives which Q must vanish and which must be
   equal or opposite; every computed value is checked against these relations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations, product

import numpy as np

from . import elements as el
from .symmetry import SymmetryResult

IRREPS = ["A1g", "Eg", "Eg", "T2g", "T2g", "T2g", "T2u", "T2u", "T2u",
          "T1u(I)", "T1u(I)", "T1u(I)", "T1u(II)", "T1u(II)", "T1u(II)"]
SITES = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]], float)


# ----------------------------------------------------------------------------
# Definitions of Q (acting on the 21-dim coordinate vector: atoms n = 0..6, components x, y, z)
# ----------------------------------------------------------------------------
def q_norm(Q) -> float:
    """|Q| = sqrt(Q2² + … + Q15²)."""
    return float(np.linalg.norm(np.asarray(Q, float)[1:15]))


DEFAULT_M0 = 1.0    # metal mass (default)
DEFAULT_M = 1.0     # donor-atom mass (default; all six donors are treated as identical)


def q_basis(m_metal: float, m_lig: float) -> np.ndarray:
    """Return the (15, 21) matrix whose row k holds the coefficients of Q_{k+1}. M0 = metal mass, M = donor mass."""
    M0, M = m_metal, m_lig
    B = np.zeros((15, 21))
    ix = {"x": 0, "y": 1, "z": 2}

    def put(row, n, c, coef):
        B[row, 3 * n + ix[c]] += coef

    s6, s12 = 1 / np.sqrt(6), 1 / np.sqrt(12)
    for n, c in ((1, "x"), (2, "y"), (3, "z")):
        put(0, n, c, s6)
    for n, c in ((4, "x"), (5, "y"), (6, "z")):
        put(0, n, c, -s6)
    # Eg
    put(1, 1, "x", .5); put(1, 2, "y", -.5); put(1, 4, "x", -.5); put(1, 5, "y", .5)
    put(2, 3, "z", 2 * s12); put(2, 6, "z", -2 * s12)
    put(2, 1, "x", -s12); put(2, 2, "y", -s12); put(2, 4, "x", s12); put(2, 5, "y", s12)
    # T2g
    put(3, 1, "y", .5); put(3, 2, "x", .5); put(3, 4, "y", -.5); put(3, 5, "x", -.5)
    put(4, 3, "x", .5); put(4, 1, "z", .5); put(4, 6, "x", -.5); put(4, 4, "z", -.5)
    put(5, 5, "z", .5); put(5, 6, "y", .5); put(5, 3, "y", -.5); put(5, 2, "z", -.5)
    # T2u
    put(6, 1, "z", .5); put(6, 2, "z", -.5); put(6, 4, "z", .5); put(6, 5, "z", -.5)
    put(7, 3, "y", .5); put(7, 6, "y", .5); put(7, 1, "y", -.5); put(7, 4, "y", -.5)
    put(8, 3, "x", .5); put(8, 6, "x", .5); put(8, 5, "x", -.5); put(8, 2, "x", -.5)
    # T1u(I)
    c1 = np.sqrt(4 * M / ((M0 + 6 * M) * (M0 + 2 * M)))
    a, b, m = c1 * (M0 + 2 * M) / (4 * np.sqrt(M)), c1 * np.sqrt(M), c1 * np.sqrt(M0)
    for row, c, four, two in ((9, "z", (1, 2, 4, 5), (3, 6)),
                              (10, "y", (1, 6, 3, 4), (2, 5)),
                              (11, "x", (5, 6, 3, 2), (1, 4))):
        for n in four:
            put(row, n, c, a)
        for n in two:
            put(row, n, c, -b)
        put(row, 0, c, -m)
    # T1u(II)
    c2 = np.sqrt(2 * M / (M0 * (M0 + 2 * M)))
    for row, c, two in ((12, "z", (3, 6)), (13, "y", (2, 5)), (14, "x", (1, 4))):
        put(row, 0, c, c2 * np.sqrt(M0))
        for n in two:
            put(row, n, c, -c2 * M0 / (2 * np.sqrt(M)))
    return B


# ----------------------------------------------------------------------------
# representative subgroups in standard orientation
# ----------------------------------------------------------------------------
def _m(rows):
    return np.array(rows, float)


E = np.eye(3)
I = -E
C4z = _m([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
C2x, C2y, C2z = np.diag([1., -1, -1]), np.diag([-1., 1, -1]), np.diag([-1., -1, 1])
S4z = np.diag([1., 1, -1]) @ C4z
C2_110 = _m([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
C2_1m10 = _m([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
SIGx, SIGz = np.diag([-1., 1, 1]), np.diag([1., 1, -1])
SIG_110 = _m([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])     # normal [110]
SIG_1m10 = _m([[0, 1, 0], [1, 0, 0], [0, 0, 1]])      # normal [1-10], i.e. the x = y plane
C3_111 = _m([[0, 0, 1], [1, 0, 0], [0, 1, 0]])


def _oh():
    out = []
    for p in permutations(range(3)):
        for s in product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i in range(3):
                M[i, p[i]] = s[i]
            out.append(M)
    return out


def _closure(gens):
    G = [E]
    frontier = [E]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                p = np.round(g @ a)
                if not any(np.array_equal(p, h) for h in G):
                    G.append(p)
                    new.append(p)
        frontier = new
    return G


GENERATORS = {
    "D4h": [C4z, C2x, I], "D2d": [S4z, C2_110], "C4v": [C4z, SIGx],
    "D2h^d": [C2z, C2_110, I], "D2h^v": [C2x, C2y, I],
    "D3d": [C3_111, C2_1m10, I], "C3v": [C3_111, SIG_1m10], "C3": [C3_111], "D3": [C3_111, C2_1m10],
    "D2": [C2z, C2_110], "C2v^v": [C2z, SIGx], "C2v^d": [C2z, SIG_110], "C2v^m": [C2_110, SIGz],
    "C2h": [C2z, I], "C2h^xy": [C2_110, I],
    "C2^z": [C2z], "C2^xy": [C2_110], "Cs^z": [SIGz], "Cs^xy": [SIG_1m10], "Ci": [I], "C1": [],
}
STANDARD_SUBGROUPS = {k: _closure(v) for k, v in GENERATORS.items()}
STANDARD_SUBGROUPS["Oh"] = _oh()


# ----------------------------------------------------------------------------
# group theory: zero / equal / opposite relations of Q implied by the standard subgroup
# ----------------------------------------------------------------------------
def _site_rep(h: np.ndarray) -> np.ndarray:
    T = np.zeros((21, 21))
    for i, p in enumerate(SITES):
        j = int(np.argmin(np.linalg.norm(SITES - h @ p, axis=1)))
        T[3 * j:3 * j + 3, 3 * i:3 * i + 3] = h
    return T


@dataclass
class Constraint:
    zero: bool                       # required to vanish by symmetry
    relation: tuple[int, int] | None = None   # (k, ±1): Q equals ±Q_k (k counted from 1)

    def text(self) -> str:
        if self.zero:
            return "0"
        if self.relation:
            k, s = self.relation
            return f"= {'' if s > 0 else '−'}Q{k}"
        return "free"


def constraints_for(label: str, m_metal=DEFAULT_M0, m_lig=DEFAULT_M) -> list[Constraint]:
    H = STANDARD_SUBGROUPS[label]
    P = sum(_site_rep(h) for h in H) / len(H)
    B = q_basis(m_metal, m_lig)
    out = []
    for k in range(15):
        zero = np.linalg.norm(P @ B[k]) < 1e-9
        rel = None
        if not zero:
            for j in range(k):
                for s in (1, -1):
                    if np.linalg.norm(P @ (B[k] - s * B[j])) < 1e-9:
                        rel = (j + 1, s)
                        break
                if rel:
                    break
        out.append(Constraint(zero, rel))
    return out


# ----------------------------------------------------------------------------
# frame selection
# ----------------------------------------------------------------------------
def _polar(X: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(X)
    return U @ Vt


def _intertwiner_projector(pairs) -> np.ndarray:
    """Orthogonal projector (acting on vec(M)) onto the 3x3 matrices M with h M = M g for all pairs.
    The dimension is that of the commutant of the standard subgroup; for approximately symmetric
    input the singular vectors of the smallest singular values are used."""
    if not pairs:
        return np.eye(9)
    hs = [h for _, h in pairs]
    A0 = np.vstack([np.kron(h, E) - np.kron(E, h.T) for h in hs])
    dim = int(np.sum(np.linalg.svd(A0, compute_uv=False) < 1e-9)) + max(0, 9 - A0.shape[0])
    A = np.vstack([np.kron(h, E) - np.kron(E, g.T) for g, h in pairs])
    _, _, Vt = np.linalg.svd(A)
    null = Vt[9 - dim:]
    return null.T @ null


def _assign(U: np.ndarray, R: np.ndarray) -> list[int]:
    """Assign the 6 donors to the 6 octahedral sites (maximum directional overlap); returns site[i] in 1..6."""
    V = U @ R.T
    best, best_s = None, -1e9
    for perm in permutations(range(6)):
        s = sum(V[i] @ SITES[perm[i] + 1] for i in range(6))
        if s > best_s:
            best, best_s = perm, s
    return [p + 1 for p in best]


def _base_frames(U: np.ndarray) -> list[np.ndarray]:
    """Initial frames from the three pseudo-trans donor pairs, combined with the 24 proper rotations of Oh."""
    ang = U @ U.T
    best, best_s = None, 9
    for perm in permutations(range(6)):
        pairs = [(perm[0], perm[1]), (perm[2], perm[3]), (perm[4], perm[5])]
        if not (perm[0] < perm[1] and perm[2] < perm[3] and perm[4] < perm[5] and perm[0] < perm[2] < perm[4]):
            continue
        s = sum(ang[a, b] for a, b in pairs)
        if s < best_s:
            best, best_s = pairs, s
    D = np.array([U[a] - U[b] for a, b in best])         # three axis directions (rows)
    R0 = _polar(D)                                        # R0 @ direction ≈ basis vector
    if np.linalg.det(R0) < 0:
        R0[2] *= -1
    return [O @ R0 for O in _oh() if np.linalg.det(O) > 0]


@dataclass
class VanVleckResult:
    label: str
    Q: np.ndarray                     # (15,) Q1..Q15 in Å
    R: np.ndarray                     # rotation: structure Cartesian frame -> standard frame
    site_of_ligand: list[int]         # site (1..6) of the i-th donor of the sphere
    ligand_at_site: list[int]         # donor index (0..5 in the sphere) at sites 1..6
    coords_std: np.ndarray            # (7,3) standard-frame coordinates: row 0 metal, row k site k
    constraints: list[Constraint]
    checks: list[str]                 # per-Q check: 'ok' / 'violated' / 'accidental-zero'
    misfit: float                     # sum |R u_i − e_site|^2 (u = unit vectors)
    zero_tol: float
    masses: tuple[float, float]
    notes: list[str] = field(default_factory=list)
    n_equivalent_frames: int = 1     # number of equally well aligned equivalent frames

    @property
    def norm(self) -> float:
        """|Q| = sqrt(Q2² + … + Q15²): length of the van Vleck vector without the A1g breathing mode Q1 (Å)."""
        return q_norm(self.Q)

    @property
    def all_ok(self) -> bool:
        return all(c != "violated" for c in self.checks)

    def to_text(self) -> str:
        rows = ["Mode\tIrrep\tQ (Å)\tConstraint\tCheck"]
        for k in range(15):
            rows.append(f"Q{k + 1}\t{IRREPS[k]}\t{self.Q[k]:.6f}\t{self.constraints[k].text()}\t{self.checks[k]}")
        rows.append(f"|Q|\t(Q2–Q15)\t{self.norm:.6f}")
        return "\n".join(rows)


def compute(coords: np.ndarray, elems: list[str], symres: SymmetryResult,
            zero_tol: float = 1e-5,
            frame_tol: float = 1e-6, ref_Q: np.ndarray | None = None,
            m_metal: float = DEFAULT_M0, m_lig: float = DEFAULT_M,
            q7_sign: int = 0, pick: int | None = None, candidates_only: bool = False):
    """coords: (7, 3) Cartesian coordinates, metal in row 0; symres: symmetry analysis of this sphere.
    ref_Q: Q of a reference structure of the same material. If given, the frame is chosen among
    the feasible, equally well aligned equivalent frames so that Q2–Q15 are closest to ref_Q.
    q7_sign: prefer frames whose Q7 has this sign. pick: use the pick-th equivalent frame.
    candidates_only: return the Q vectors of all equivalent frames instead of a result."""
    X = np.asarray(coords, float)
    rel = X[1:] - X[0]
    U = rel / np.linalg.norm(rel, axis=1)[:, None]
    label = symres.label
    notes = []
    if label not in STANDARD_SUBGROUPS:
        raise ValueError(f"No standard orientation defined for {label}")
    H = STANDARD_SUBGROUPS[label]
    G = [o.matrix for o in symres.ops]

    lig_elems = {e for e in elems[1:]}
    m_met = float(m_metal)          # user-defined masses; all six donors treated as identical
    m_lig = float(m_lig)
    QB = q_basis(m_met, m_lig)

    def q_of(R, site):
        cs = np.zeros((7, 3))
        for i in range(6):
            cs[site[i]] = R @ rel[i]
        return QB @ cs.reshape(21)

    cands = []          # all feasible frames: (misfit, R, site)
    for R in _base_frames(U):
        for _ in range(12):
            # 1) pair the molecule's operations with elements of the standard subgroup
            pairs, used, ok = [], set(), True
            for g in G:
                Rg = R @ g @ R.T
                d = [np.linalg.norm(Rg - h) for h in H]
                j = int(np.argmin(d))
                if d[j] > 0.9 or j in used:
                    ok = False
                    break
                used.add(j)
                pairs.append((g, H[j]))
            if not ok or len(used) != len(H):
                break
            # 2) assign donors to sites, 3) best rotation under the constraint
            site = _assign(U, R)
            Bm = sum(np.outer(SITES[site[i]], U[i]) for i in range(6))   # max tr(Rᵀ Bm)
            Pm = _intertwiner_projector(pairs)
            Rn = _polar((Pm @ Bm.reshape(9)).reshape(3, 3))
            if np.linalg.det(Rn) < 0:
                Rn = -Rn
            if np.linalg.norm(Rn - R) < 1e-12:
                R = Rn
                break
            R = Rn
        if not ok:
            continue
        # final check: symmetry elements exactly in standard positions
        if any(min(np.linalg.norm(R @ g @ R.T - h) for h in H) > frame_tol for g in G):
            continue
        site = _assign(U, R)
        misfit = float(sum(np.sum((R @ U[i] - SITES[site[i]]) ** 2) for i in range(6)))
        if not any(np.linalg.norm(R - c[1]) < 1e-6 for c in cands):
            cands.append((misfit, R, site))
    if not cands:
        raise RuntimeError("Could not find a coordinate frame placing the symmetry elements in standard positions.")
    mmin = min(c[0] for c in cands)
    good = [c for c in cands if c[0] <= mmin + 1e-6]          # equivalent frames with the best octahedral alignment
    good.sort(key=lambda c: float(np.linalg.norm(c[1] - E)))  # fixed order, so that `pick` is reproducible
    if candidates_only:                                       # only return Q for each equivalent frame
        return [q_of(c[1], c[2]) for c in good]
    # q7_sign = ±1: prefer equivalent frames whose Q7 has that sign (no effect if Q7 = 0)
    def sign_penalty(c):
        if not q7_sign:
            return 0
        q7 = q_of(c[1], c[2])[6]
        return 1 if (abs(q7) >= zero_tol and np.sign(q7) != q7_sign) else 0

    if ref_Q is not None:
        ref_Q = np.asarray(ref_Q, float)
        key = lambda c: (sign_penalty(c), float(np.linalg.norm(q_of(c[1], c[2])[1:] - ref_Q[1:])),
                         float(np.linalg.norm(c[1] - E)))
    else:
        key = lambda c: (sign_penalty(c), float(np.linalg.norm(c[1] - E)))
    misfit, R, site = good[pick] if pick is not None else min(good, key=key)
    n_equiv = len(good)

    coords_std = np.zeros((7, 3))
    ligand_at_site = [0] * 6
    for i in range(6):
        coords_std[site[i]] = R @ rel[i]
        ligand_at_site[site[i] - 1] = i

    if len(lig_elems) > 1:
        notes.append(f"Mixed donor atoms ({', '.join(sorted(lig_elems))}); all six are treated as identical "
                     f"atoms of mass M = {m_lig:g} for Q.")
    Q = QB @ coords_std.reshape(21)

    cons = constraints_for(label, m_met, m_lig)
    checks = []
    for k, c in enumerate(cons):
        if c.zero:
            checks.append("ok" if abs(Q[k]) < zero_tol else "violated")
        elif c.relation:
            j, s = c.relation
            checks.append("ok" if abs(Q[k] - s * Q[j - 1]) < zero_tol else "violated")
        else:
            checks.append("accidental-zero" if (abs(Q[k]) < zero_tol and k > 0) else "ok")
    return VanVleckResult(label, Q, R, site, ligand_at_site, coords_std, cons, checks, misfit,
                          zero_tol, (m_met, m_lig), notes, n_equiv)

