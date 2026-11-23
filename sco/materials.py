"""Decide whether two complexes are the "same material" (same chemical structure) for comparing Q.

Criterion: the metal-containing molecules have the same non-H composition (summed with site
occupancies, so disorder is handled) and the same donor elements.
"""
from __future__ import annotations

from collections import defaultdict

from .coordination import CoordinationSphere


def molecule_key(sphere: CoordinationSphere) -> tuple:
    st = sphere.structure
    mol = st.molecule_of[sphere.metal_index]
    comp = defaultdict(float)
    for i, a in enumerate(st.atoms):
        if st.molecule_of[i] == mol and a.element != "H":
            comp[a.element] += a.occupancy
    donors = tuple(sorted(a.element for a in sphere.ligands))
    return tuple(sorted((k, round(v, 1)) for k, v in comp.items())), sphere.metal.element, donors


def formula_text(key: tuple) -> str:
    comp, metal, donors = key
    f = " ".join(f"{e}{v:g}" if v != 1 else e for e, v in comp)
    return f"{f}  [{metal}{''.join(donors)}]"


def same_material(k1: tuple, k2: tuple, tol: float = 0.3) -> bool:
    c1, m1, d1 = k1
    c2, m2, d2 = k2
    if m1 != m2 or d1 != d2:
        return False
    a, b = dict(c1), dict(c2)
    if set(a) != set(b):
        return False
    return all(abs(a[e] - b[e]) <= tol for e in a)
