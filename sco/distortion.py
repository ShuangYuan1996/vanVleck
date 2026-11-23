"""Octahedral distortion parameters, following the definitions used by OctaDist
(R. Ketkaew et al., Dalton Trans., 2021, 50, 1086). This is an independent implementation.

<D>  = d_mean                                    mean M–L distance (Å)
zeta = sum_{i=1..6} |d_i − d_mean|               bond-length distortion (Å)
Delta = (1/6) sum_{i=1..6} ((d_i − d_mean)/d_mean)^2   (dimensionless; "tilting distortion" in OctaDist)
Sigma = sum_{i=1..12} |phi_i − 90°|              12 cis L–M–L angles (°)
Theta = sum_{i=1..24} |theta_i − 60°|            24 torsion angles over the 4 pairs of opposite faces (°)

Theta follows OctaDist: the 8 faces form 4 pairs of opposite faces. For each pair, one face is the
reference; the metal and the 3 donors of the opposite face are projected orthogonally onto the
reference plane, and the 6 signed angles between consecutive "reference donor -> opposite donor"
vectors are taken counter-clockwise (ideal value 60°). The same is done with the other face as
reference and the two results are averaged; the 4 pairs are summed.
Equivalently, with t_1…t_8 the sums obtained with each face as reference, Theta = sum(t)/2,
Theta_min = sum of the 4 smallest t, Theta_max = sum of the 4 largest t (as in OctaDist).
V = volume of the convex hull (octahedron) of the 6 donor atoms (Å³).
rmsd() gives the RMSD between two coordination spheres after Kabsch superposition.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations

import numpy as np


@dataclass
class Distortion:
    zeta: float          # Å
    sigma: float         # °
    theta: float         # °
    d_mean: float        # Å
    delta: float = 0.0   # dimensionless
    theta_min: float = 0.0   # °
    theta_max: float = 0.0   # °
    volume: float = 0.0  # Å³


def _trans_pairs(U: np.ndarray) -> list[tuple[int, int]]:
    """Split the 6 donors into 3 pseudo-trans pairs (maximising the three trans angles)."""
    G = U @ U.T
    best, best_s = None, 9.0
    for p in permutations(range(6)):
        if not (p[0] < p[1] and p[2] < p[3] and p[4] < p[5] and p[0] < p[2] < p[4]):
            continue
        s = G[p[0], p[1]] + G[p[2], p[3]] + G[p[4], p[5]]
        if s < best_s:
            best, best_s = [(p[0], p[1]), (p[2], p[3]), (p[4], p[5])], s
    return best


def _signed_angle(a: np.ndarray, b: np.ndarray, n: np.ndarray) -> float:
    """Signed in-plane angle from a to b (°), counter-clockwise about the normal n."""
    return float(np.degrees(np.arctan2(np.dot(np.cross(a, b), n), np.dot(a, b))))


def _face_twist(L: np.ndarray, ref: tuple, opp_of: dict) -> float:
    """With face ref (three donor indices) as reference, return the sum of the 6 terms |theta − 60|."""
    P = L[list(ref)]
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n /= np.linalg.norm(n)
    c = P.mean(axis=0)
    if np.dot(n, c) < 0:            # normal points away from the metal (metal at the origin)
        n = -n
    proj = lambda x: x - np.dot(x - c, n) * n
    m = proj(np.zeros(3))           # projected metal position
    vec = {i: proj(L[i]) - m for i in range(6)}
    # order the reference donors counter-clockwise
    e1 = vec[ref[0]] / np.linalg.norm(vec[ref[0]])
    order = sorted(ref, key=lambda i: _signed_angle(e1, vec[i], n) % 360.0)
    a, b, cc = order
    # ideally the donor trans to cc lies between a and b, and so on
    seq = [a, opp_of[cc], b, opp_of[a], cc, opp_of[b]]
    total = 0.0
    for k in range(6):
        th = _signed_angle(vec[seq[k]], vec[seq[(k + 1) % 6]], n)
        total += abs(th - 60.0)
    return total


def compute(coords: np.ndarray) -> Distortion:
    """coords: (7, 3); row 0 is the metal, rows 1–6 the donors (any order)."""
    X = np.asarray(coords, float)
    L = X[1:] - X[0]
    d = np.linalg.norm(L, axis=1)
    d_mean = float(d.mean())
    zeta = float(np.sum(np.abs(d - d_mean)))

    U = L / d[:, None]
    ang = sorted(float(np.degrees(np.arccos(np.clip(U[i] @ U[j], -1, 1))))
                 for i, j in combinations(range(6), 2))
    sigma = float(sum(abs(a - 90.0) for a in ang[:12]))      # drop the 3 largest (trans) angles

    delta = float(np.mean(((d - d_mean) / d_mean) ** 2))

    pairs = _trans_pairs(U)
    opp_of = {}
    for i, j in pairs:
        opp_of[i], opp_of[j] = j, i
    # the 8 faces: one donor from each trans pair
    eight = [_face_twist(L, tuple(pairs[k][(pick >> k) & 1] for k in range(3)), opp_of) for pick in range(8)]
    eight.sort()
    theta = float(sum(eight) / 2.0)
    theta_min = float(sum(eight[:4]))
    theta_max = float(sum(eight[4:]))
    return Distortion(zeta, sigma, theta, d_mean, delta, theta_min, theta_max, hull_volume(L))


def hull_volume(P: np.ndarray) -> float:
    """Convex-hull volume of a few points: sum of tetrahedra (centroid + hull face) over all hull faces."""
    P = np.asarray(P, float)
    c = P.mean(axis=0)
    vol = 0.0
    n = len(P)
    for i, j, k in combinations(range(n), 3):
        nrm = np.cross(P[j] - P[i], P[k] - P[i])
        if np.linalg.norm(nrm) < 1e-12:
            continue
        side = [(P[m] - P[i]) @ nrm for m in range(n) if m not in (i, j, k)]
        if all(x <= 1e-9 for x in side) or all(x >= -1e-9 for x in side):
            vol += abs(np.dot(P[i] - c, nrm)) / 6.0
    return float(vol)


def rmsd(A: np.ndarray, B: np.ndarray) -> float:
    """RMSD (Å) between two ordered coordinate sets after centring and optimal (Kabsch) rotation."""
    A = np.asarray(A, float) - np.asarray(A, float).mean(axis=0)
    B = np.asarray(B, float) - np.asarray(B, float).mean(axis=0)
    H = A.T @ B
    U_, S, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U_.T))])
    R = Vt.T @ D @ U_.T
    return float(np.sqrt(np.mean(np.sum((A @ R.T - B) ** 2, axis=1))))


def sci(x: float, digits: int = 3) -> str:
    """Readable scientific notation, e.g. 6.07×10⁻⁴."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m:.{digits - 1}f}×10{sup}"

