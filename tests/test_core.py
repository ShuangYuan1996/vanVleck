"""Regression tests for vanVleck.

Run with   python -m pytest tests      (or simply   python tests/test_core.py)
Only numpy is needed; the GUI is not imported.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from sco import distortion as dist                      # noqa: E402
from sco.structure import load_file, parse_symop, format_symop   # noqa: E402
from sco.coordination import find_coordination_sphere    # noqa: E402
from sco.symmetry import analyze                          # noqa: E402
from sco.vanvleck import (STANDARD_SUBGROUPS, compute, constraints_for, q_basis,   # noqa: E402
                          _site_rep)

EXAMPLE = os.path.join(ROOT, "examples", "synthetic_FeN6.xyz")
IDEAL = np.vstack([[0, 0, 0], np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                                        [-1, 0, 0], [0, -1, 0], [0, 0, -1.]]) * 2.0])
ELEMS = ["Fe"] + ["N"] * 6


def _random_rotation(rng):
    q = rng.normal(size=4)
    a, b, c, d = q / np.linalg.norm(q)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]])


def test_symop_roundtrip():
    R, t = parse_symop("-x+1/2, y, -z+1/2")
    assert np.allclose(R, np.diag([-1, 1, -1])) and np.allclose(t, [0.5, 0, 0.5])
    assert format_symop(R, t) == "1/2-x,y,1/2-z"


def test_ideal_octahedron():
    sym = analyze(IDEAL, ELEMS, 1e-5)
    assert sym.label == "Oh" and sym.order == 48
    r = compute(IDEAL, ELEMS, sym)
    assert np.isclose(r.Q[0], 3 * 2.0 * 2 / np.sqrt(6))        # Q1 = 6 d / sqrt(6)
    assert np.allclose(r.Q[1:], 0, atol=1e-12) and r.norm < 1e-12
    d = dist.compute(IDEAL)
    assert d.zeta < 1e-12 and d.delta < 1e-12 and d.sigma < 1e-9 and d.theta < 1e-9
    assert np.isclose(d.volume, 4 / 3 * 2.0 ** 3)


def test_trigonal_prism_theta():
    h, r = 1.2, 1.6
    pts = [[0, 0, 0]] + [[r * np.cos(a), r * np.sin(a), s * h]
                         for s in (1, -1) for a in np.radians([0, 120, 240])]
    assert dist.compute(np.array(pts)).theta > 1000            # close to the trigonal-prism maximum


def test_q_basis_orthonormal_equal_masses():
    B = q_basis(1.0, 1.0)
    assert np.allclose(B @ B.T, np.eye(15), atol=1e-12)


def test_synthetic_examples():
    """Every frame of the example file has the point group named in its comment line,
    and all zero / equal / opposite relations of Q hold."""
    structures = load_file(EXAMPLE)
    assert len(structures) == 19
    for st in structures:
        expected = st.name.split()[-1]
        sp = find_coordination_sphere(st, st.metal_indices()[0])
        elems = [sp.metal.element] + [a.element for a in sp.ligands]
        sym = analyze(sp.coordinates(), elems, 1e-5)
        assert sym.label == expected, (expected, sym.label)
        r = compute(sp.coordinates(), elems, sym, zero_tol=1e-5)
        assert r.all_ok, (expected, r.checks)


def test_all_subgroups_random_orientation():
    """Random distortions of each standard subgroup, in random orientations, are recognised,
    placed in the standard frame and satisfy every symmetry relation."""
    rng = np.random.default_rng(7)
    for label, H in STANDARD_SUBGROUPS.items():
        P = sum(_site_rep(h) for h in H) / len(H)
        for _ in range(2):
            cfg = IDEAL + 0.12 * (P @ rng.normal(size=21)).reshape(7, 3)
            cfg -= cfg[0]
            X = cfg @ _random_rotation(rng).T + rng.normal(size=3) * 3
            sym = analyze(X, ELEMS, 1e-5)
            assert sym.label == label, (label, sym.label)
            r = compute(X, ELEMS, sym, zero_tol=1e-5)
            assert r.all_ok, (label, r.checks)
            for k, c in enumerate(constraints_for(label)):
                if c.zero:
                    assert abs(r.Q[k]) < 1e-6


def test_rmsd_rotation_invariant():
    rng = np.random.default_rng(1)
    A = IDEAL + 0.1 * rng.normal(size=(7, 3))
    B = A @ _random_rotation(rng).T + 2.0
    assert dist.rmsd(A, B) < 1e-10


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("passed", name)
