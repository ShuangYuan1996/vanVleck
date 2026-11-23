"""Inner coordination sphere: the metal centre plus its 6 nearest donor atoms."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .structure import Atom, Structure

MAX_DONOR_DISTANCE = 3.0   # donors farther than this are flagged as suspicious


@dataclass
class CoordinationSphere:
    structure: Structure
    metal_index: int
    ligand_indices: list[int]          # sorted by increasing M–L distance
    distances: list[float]

    @property
    def metal(self) -> Atom:
        return self.structure.atoms[self.metal_index]

    @property
    def ligands(self) -> list[Atom]:
        return [self.structure.atoms[i] for i in self.ligand_indices]

    @property
    def indices(self) -> list[int]:
        return [self.metal_index] + list(self.ligand_indices)

    @property
    def metal_coord(self) -> np.ndarray:
        return self.metal.cart.copy()

    @property
    def ligand_coords(self) -> np.ndarray:
        """(6, 3) Cartesian coordinates of the donor atoms (Å)."""
        return np.array([a.cart for a in self.ligands])

    def coordinates(self) -> np.ndarray:
        """(7, 3): row 0 is the metal, rows 1–6 the donor atoms."""
        return np.vstack([self.metal_coord, self.ligand_coords])

    @property
    def mean_distance(self) -> float:
        return float(np.mean(self.distances))

    @property
    def suspicious(self) -> bool:
        return len(self.ligand_indices) < 6 or max(self.distances) > MAX_DONOR_DISTANCE

    def cis_faces(self, max_cis_angle: float = 135.0) -> list[tuple[int, int, int]]:
        """The 8 faces of a pseudo-octahedron: triples of mutually cis donors (L–M–L < max_cis_angle)."""
        m = self.metal_coord
        v = self.ligand_coords - m
        v /= np.linalg.norm(v, axis=1)[:, None]
        ang = np.degrees(np.arccos(np.clip(v @ v.T, -1, 1)))
        n = len(v)
        faces = []
        for a in range(n):
            for b in range(a + 1, n):
                for c in range(b + 1, n):
                    if ang[a, b] < max_cis_angle and ang[a, c] < max_cis_angle and ang[b, c] < max_cis_angle:
                        li = self.ligand_indices
                        faces.append((li[a], li[b], li[c]))
        return faces

    def to_text(self) -> str:
        """Tab-separated text, ready to paste into Excel / Origin."""
        st = self.structure
        rows = ["Atom\tElement\tx(Å)\ty(Å)\tz(Å)\tfrac_x\tfrac_y\tfrac_z\td(M-L)(Å)\tSymmetry"]
        for idx, d in [(self.metal_index, 0.0)] + list(zip(self.ligand_indices, self.distances)):
            a = st.atoms[idx]
            f = a.frac if a.frac is not None else [np.nan] * 3
            rows.append(f"{a.label}\t{a.element}\t{a.cart[0]:.5f}\t{a.cart[1]:.5f}\t{a.cart[2]:.5f}\t"
                        f"{f[0]:.5f}\t{f[1]:.5f}\t{f[2]:.5f}\t"
                        f"{'' if idx == self.metal_index else f'{d:.4f}'}\t{a.symop}")
        return "\n".join(rows)


def find_coordination_sphere(st: Structure, metal_index: int, n: int = 6) -> CoordinationSphere:
    """Take the atoms bonded to the metal first; if fewer than n, fill up with the nearest non-H atoms."""
    xyz = st.coords
    m = xyz[metal_index]
    d_all = np.linalg.norm(xyz - m, axis=1)
    bonded = [j for j in st.neighbors(metal_index) if st.atoms[j].element != "H"]
    bonded.sort(key=lambda j: d_all[j])
    chosen = bonded[:n]
    if len(chosen) < n:
        order = np.argsort(d_all)
        for j in order:
            j = int(j)
            if j == metal_index or j in chosen or st.atoms[j].element == "H":
                continue
            chosen.append(j)
            if len(chosen) == n:
                break
    chosen.sort(key=lambda j: d_all[j])
    return CoordinationSphere(st, metal_index, chosen, [float(d_all[j]) for j in chosen])


def sphere_for_clicked_atom(st: Structure, atom_index: int) -> CoordinationSphere | None:
    """For a clicked atom, return the coordination sphere of the nearest metal in the same molecule."""
    metals = st.metal_indices()
    if not metals:
        return None
    mol = st.molecule_of[atom_index] if st.molecule_of else 0
    same = [m for m in metals if st.molecule_of[m] == mol]
    if not same:
        return None
    xyz = st.coords
    best = min(same, key=lambda m: np.linalg.norm(xyz[m] - xyz[atom_index]))
    return find_coordination_sphere(st, best)
