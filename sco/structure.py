"""Crystal/molecular structure model: CIF and XYZ reading, symmetry expansion, bonding, molecule completion."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

from . import elements as el
from .cif_parser import CifBlock, cif_float, parse_cif

BOND_TOLERANCE = 0.40      # bond if d < r1 + r2 + 0.40 Å (as in the CSD / Mercury)
MAX_GROWN_ATOMS = 6000     # growth limit for polymeric (infinite) structures


# ----------------------------------------------------------------------------
# data classes
# ----------------------------------------------------------------------------
@dataclass
class Atom:
    label: str
    element: str
    cart: np.ndarray                 # Cartesian coordinates (Å)
    frac: np.ndarray | None = None   # fractional coordinates (CIF only)
    occupancy: float = 1.0
    disorder_group: str = "."
    symop: str = "x,y,z"             # symmetry operator (incl. translation) that generated the atom
    asym_index: int = -1             # index of the atom in the asymmetric unit


@dataclass
class Cell:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    @property
    def matrix(self) -> np.ndarray:
        """Columns are the Cartesian components of a, b, c (a along x, b in the xy plane)."""
        al, be, ga = np.radians([self.alpha, self.beta, self.gamma])
        ca, cb, cg, sg = np.cos(al), np.cos(be), np.cos(ga), np.sin(ga)
        cy = (ca - cb * cg) / sg
        cz = np.sqrt(max(1.0 - cb * cb - cy * cy, 0.0))
        return np.array([[self.a, self.b * cg, self.c * cb],
                         [0.0, self.b * sg, self.c * cy],
                         [0.0, 0.0, self.c * cz]])

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.matrix)))


@dataclass
class Structure:
    name: str
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[tuple[int, int]] = field(default_factory=list)
    molecule_of: list[int] = field(default_factory=list)   # molecule index of each atom
    cell: Cell | None = None
    space_group: str = ""
    symops: list[str] = field(default_factory=list)
    temperature: float | None = None
    formula: str = ""
    compound_name: str = ""
    ccdc: str = ""
    formula_sum: str = ""
    source_file: str = ""
    source_format: str = ""
    notes: list[str] = field(default_factory=list)

    # ---- convenience properties ----
    @property
    def coords(self) -> np.ndarray:
        if not self.atoms:
            return np.zeros((0, 3))
        return np.array([a.cart for a in self.atoms])

    @property
    def n_molecules(self) -> int:
        return (max(self.molecule_of) + 1) if self.molecule_of else 0

    def neighbors(self, idx: int) -> list[int]:
        return [j if i == idx else i for i, j in self.bonds if idx in (i, j)]

    def metal_indices(self) -> list[int]:
        tm = [i for i, a in enumerate(self.atoms) if a.element in el.TRANSITION_METALS]
        if tm:
            return tm
        return [i for i, a in enumerate(self.atoms) if a.element in el.METALS]

    @property
    def short_name(self) -> str:
        """For path-like block names such as 'c:\\datasets\\leeds\\je40a1~1' keep only the last part."""
        n = re.split(r"[\\/]", self.name)[-1] or self.name
        return n if len(n) <= 24 else "…" + n[-23:]   # keep the end (often holds the temperature)

    @property
    def short_compound(self) -> str:
        """'... diperchlorate (compound 5)' -> 'compound 5'; otherwise the first 40 characters."""
        m = re.search(r"\((compound\s+[^)]+)\)", self.compound_name, re.I)
        if m:
            return m.group(1)
        if self.compound_name:
            return self.compound_name[:40] + ("…" if len(self.compound_name) > 40 else "")
        return self.formula

    def display_name(self) -> str:
        parts = [p for p in (self.space_group,
                             f"{self.temperature:g} K" if self.temperature is not None else "") if p]
        txt = self.short_name
        if self.short_compound:
            txt += f" — {self.short_compound}"
        return txt + (f"  ({', '.join(parts)})" if parts else "")


# ----------------------------------------------------------------------------
# symmetry operators
# ----------------------------------------------------------------------------
def parse_symop(text: str) -> tuple[np.ndarray, np.ndarray]:
    """'-x+1/2, y, -z+1/2' -> (R 3x3, t 3)."""
    s = text.replace(" ", "").replace("'", "").lower()
    parts = s.split(",")
    if len(parts) != 3:
        raise ValueError(f"Cannot parse symmetry operator: {text}")
    R = np.zeros((3, 3))
    t = np.zeros(3)
    for row, expr in enumerate(parts):
        for sign, term in re.findall(r"([+-]?)([^+-]+)", expr):
            sg = -1.0 if sign == "-" else 1.0
            term = term.replace("*", "")
            if term[-1] in "xyz":
                coef = term[:-1]
                c = float(Fraction(coef)) if coef else 1.0
                R[row, "xyz".index(term[-1])] += sg * c
            else:
                t[row] += sg * float(Fraction(term))
    return R, t


def format_symop(R: np.ndarray, t: np.ndarray) -> str:
    """Format (R, t) as e.g. '1-x,-y,1/2+z'."""
    out = []
    for row in range(3):
        s = ""
        tv = Fraction(float(t[row])).limit_denominator(12)
        if tv != 0:
            s += str(tv)
        for col, v in enumerate("xyz"):
            c = int(round(R[row, col]))
            if c == 0:
                continue
            if c == 1:
                s += ("+" if s else "") + v
            elif c == -1:
                s += "-" + v
            else:
                s += (f"+{c}" if (s and c > 0) else str(c)) + v
        out.append(s or "0")
    return ",".join(out)


# ----------------------------------------------------------------------------
# bonding
# ----------------------------------------------------------------------------
def _disorder_compatible(g1: str, g2: str) -> bool:
    return g1 in (".", "0", "?") or g2 in (".", "0", "?") or g1 == g2


def _bond_cutoff(e1: str, e2: str) -> float:
    return el.covalent_radius(e1) + el.covalent_radius(e2) + BOND_TOLERANCE


def compute_bonds(atoms: list[Atom]) -> list[tuple[int, int]]:
    """Distance-based bonding for non-periodic systems."""
    n = len(atoms)
    if n < 2:
        return []
    xyz = np.array([a.cart for a in atoms])
    rad = np.array([el.covalent_radius(a.element) for a in atoms])
    is_h = np.array([a.element == "H" for a in atoms])
    bonds = []
    for i in range(n - 1):
        d = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1)
        cut = rad[i] + rad[i + 1:] + BOND_TOLERANCE
        ok = (d < cut) & (d > 0.4)
        if is_h[i]:
            ok &= ~is_h[i + 1:]
        for k in np.nonzero(ok)[0]:
            j = i + 1 + int(k)
            if _disorder_compatible(atoms[i].disorder_group, atoms[j].disorder_group):
                bonds.append((i, j))
    return bonds


def assign_molecules(n_atoms: int, bonds: list[tuple[int, int]]) -> list[int]:
    parent = list(range(n_atoms))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in bonds:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    ids, out = {}, []
    for i in range(n_atoms):
        r = find(i)
        out.append(ids.setdefault(r, len(ids)))
    return out


# ----------------------------------------------------------------------------
# CIF -> Structure
# ----------------------------------------------------------------------------
def _cell_from_block(b: CifBlock) -> Cell | None:
    vals = [cif_float(b.get(t)) for t in (
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")]
    if any(v is None for v in vals[:3]):
        return None
    vals = [v if v is not None else 90.0 for v in vals]
    return Cell(*vals)


def _symops_from_block(b: CifBlock) -> list[str]:
    tags = ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
    loop = b.find_loop(*tags)
    if loop:
        ops = b.column(loop, *tags)
        if ops:
            return list(ops)
    return ["x,y,z"]


def structure_from_block(b: CifBlock, source_file: str = "") -> Structure | None:
    loop = b.find_loop("_atom_site_fract_x")
    cell = _cell_from_block(b)
    if loop is None or cell is None:
        return None
    col = lambda *t: b.column(loop, *t)
    labels = col("_atom_site_label")
    fx, fy, fz = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    types = col("_atom_site_type_symbol")
    occs = col("_atom_site_occupancy")
    groups = col("_atom_site_disorder_group")
    n = len(fx)
    if labels is None:
        labels = [f"X{i + 1}" for i in range(n)]

    asym = []   # (label, element, frac, occ, group)
    for i in range(n):
        f = [cif_float(fx[i]), cif_float(fy[i]), cif_float(fz[i])]
        if any(v is None for v in f):
            continue
        sym = el.normalize_symbol(types[i]) if types else el.normalize_symbol(labels[i])
        if not el.is_element(sym):
            sym = el.normalize_symbol(labels[i])
        occ = cif_float(occs[i]) if occs else None
        grp = groups[i] if groups else "."
        asym.append((labels[i], sym, np.array(f), 1.0 if occ is None else occ, grp))
    if not asym:
        return None

    st = Structure(name=b.name, cell=cell, source_file=source_file, source_format="CIF")
    st.space_group = (b.get("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m",
                            default="") or "").strip()
    st.symops = _symops_from_block(b)
    st.temperature = cif_float(b.get("_diffrn_ambient_temperature", "_cell_measurement_temperature"))
    st.formula = (b.get("_chemical_formula_moiety", "_chemical_formula_sum", default="") or "").strip()
    st.formula_sum = (b.get("_chemical_formula_sum", default="") or "").strip()
    st.ccdc = (b.get("_database_code_depnum_ccdc_archive", default="") or "").replace("CCDC", "").strip()
    st.compound_name = " ".join((b.get("_chemical_name_systematic", "_chemical_name_common",
                                       default="") or "").split())
    _build_molecules(st, asym)
    return st


def _build_molecules(st: Structure, asym: list) -> None:
    """Expand the asymmetric unit with the symmetry operators and grow complete molecules along bonds (like "complete molecules" in Mercury)."""
    M = st.cell.matrix
    ops = []
    for s in st.symops:
        try:
            ops.append(parse_symop(s))
        except Exception:
            st.notes.append(f"Ignored unparsable symmetry operator: {s}")
    id_idx = next((k for k, (R, t) in enumerate(ops)
                   if np.allclose(R, np.eye(3)) and np.allclose(t % 1.0, 0)), None)
    if id_idx is None:
        ops.insert(0, (np.eye(3), np.zeros(3)))
        id_idx = 0

    # 1) all atoms in the unit cell (P1 expansion), with their origin
    p1_frac, p1_src = [], []      # src = (asym_idx, op_idx, wrapping translation n0)
    identity_p1 = {}
    for ai, (_, _, f, _, _) in enumerate(asym):
        mine = []
        for k, (R, t) in enumerate(ops):
            g = R @ f + t
            n0 = -np.floor(g + 1e-8)
            w = g + n0
            dup = False
            for pi in mine:
                d = w - p1_frac[pi]
                d -= np.round(d)
                if np.linalg.norm(M @ d) < 0.2:
                    dup = True
                    if k == id_idx:
                        identity_p1[ai] = pi
                    break
            if dup:
                continue
            mine.append(len(p1_frac))
            if k == id_idx:
                identity_p1[ai] = len(p1_frac)
            p1_frac.append(w)
            p1_src.append((ai, k, n0))
    p1_frac = np.array(p1_frac)
    N = len(p1_frac)
    elem = [asym[s[0]][1] for s in p1_src]
    grp = [asym[s[0]][4] for s in p1_src]
    rad = np.array([el.covalent_radius(e) for e in elem])
    is_h = np.array([e == "H" for e in elem])

    # 2) periodic neighbour list: nbr[i] = [(j, lattice translation)]
    shifts = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)], float)
    nbr = [[] for _ in range(N)]
    for i in range(N):
        diff = p1_frac[None, :, :] - p1_frac[i][None, None, :] + shifts[:, None, :]  # (27,N,3)
        d = np.linalg.norm(diff @ M.T, axis=2)
        cut = rad[i] + rad[None, :] + BOND_TOLERANCE
        ok = (d < cut) & (d > 0.4)
        if is_h[i]:
            ok &= ~is_h[None, :]
        for s_idx, j in zip(*np.nonzero(ok)):
            if _disorder_compatible(grp[i], grp[j]):
                nbr[i].append((int(j), tuple(int(v) for v in shifts[s_idx])))

    # 3) breadth-first growth starting from the asymmetric unit
    nodes, node_index = [], {}
    covered_asym = set()
    truncated = False
    for ai in range(len(asym)):
        if ai in covered_asym:
            continue
        pi = identity_p1[ai]
        n_start = tuple(int(v) for v in np.round(asym[ai][2] - p1_frac[pi]))
        start = (pi, n_start)
        if start in node_index:
            continue
        queue = [start]
        node_index[start] = len(nodes)
        nodes.append(start)
        while queue:
            i, n = queue.pop()
            covered_asym.add(p1_src[i][0])
            for j, s in nbr[i]:
                key = (j, (n[0] + s[0], n[1] + s[1], n[2] + s[2]))
                if key in node_index:
                    continue
                if len(nodes) >= MAX_GROWN_ATOMS:
                    truncated = True
                    continue
                node_index[key] = len(nodes)
                nodes.append(key)
                queue.append(key)
    if truncated:
        st.notes.append(f"Structure appears to be polymeric; growth truncated at {MAX_GROWN_ATOMS} atoms.")

    # 4) build atoms and bonds
    for i, n in nodes:
        ai, k, n0 = p1_src[i]
        label, e, _, occ, g = asym[ai]
        frac = p1_frac[i] + np.array(n, float)
        R, t = ops[k]
        st.atoms.append(Atom(label=label, element=e, cart=M @ frac, frac=frac,
                             occupancy=occ, disorder_group=g,
                             symop=format_symop(R, t + n0 + np.array(n, float)),
                             asym_index=ai))
    bonds = set()
    for (i, n), a in node_index.items():
        for j, s in nbr[i]:
            b = node_index.get((j, (n[0] + s[0], n[1] + s[1], n[2] + s[2])))
            if b is not None and a != b:
                bonds.add((min(a, b), max(a, b)))
    st.bonds = sorted(bonds)
    st.molecule_of = assign_molecules(len(st.atoms), st.bonds)


# ----------------------------------------------------------------------------
# XYZ
# ----------------------------------------------------------------------------
def read_xyz(text: str, source_file: str = "") -> list[Structure]:
    lines = text.splitlines()
    out, i, frame = [], 0, 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            raise ValueError(f"XYZ line {i + 1}: expected the number of atoms")
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""
        atoms, counts = [], {}
        for k in range(n):
            parts = lines[i + 2 + k].split()
            sym = el.normalize_symbol(parts[0])
            counts[sym] = counts.get(sym, 0) + 1
            atoms.append(Atom(label=f"{sym}{counts[sym]}", element=sym,
                              cart=np.array([float(parts[1]), float(parts[2]), float(parts[3])])))
        frame += 1
        name = f"Frame {frame}" + (f": {comment[:50]}" if comment else "")
        st = Structure(name=name, atoms=atoms, source_file=source_file, source_format="XYZ")
        st.formula = " ".join(f"{e}{c}" for e, c in sorted(counts.items()))
        st.bonds = compute_bonds(atoms)
        st.molecule_of = assign_molecules(len(atoms), st.bonds)
        out.append(st)
        i += 2 + n
    return out


# ----------------------------------------------------------------------------
# entry point
# ----------------------------------------------------------------------------
def _read_text(path: str) -> str:
    with open(path, "rb") as fh:
        raw = fh.read()
    for enc in ("utf-8-sig", "gbk", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def load_file(path: str) -> list[Structure]:
    """Read a .cif or .xyz file and return every dataset that contains atomic coordinates."""
    text = _read_text(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xyz":
        return read_xyz(text, path)
    structures = []
    for blk in parse_cif(text):
        st = structure_from_block(blk, path)
        if st is not None:
            structures.append(st)
    return structures
