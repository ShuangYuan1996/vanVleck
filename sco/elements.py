"""Element data: colours (Jmol/CPK), CSD covalent radii and van der Waals radii."""

# symbol: (colour, covalent radius Å, van der Waals radius Å)
_DATA = {
    "H": ("#FFFFFF", 0.23, 1.10), "He": ("#D9FFFF", 1.50, 1.40),
    "Li": ("#CC80FF", 1.28, 1.82), "Be": ("#C2FF00", 0.96, 1.53),
    "B": ("#FFB5B5", 0.83, 1.92), "C": ("#909090", 0.68, 1.70),
    "N": ("#3050F8", 0.68, 1.55), "O": ("#FF0D0D", 0.68, 1.52),
    "F": ("#90E050", 0.64, 1.47), "Ne": ("#B3E3F5", 1.50, 1.54),
    "Na": ("#AB5CF2", 1.66, 2.27), "Mg": ("#8AFF00", 1.41, 1.73),
    "Al": ("#BFA6A6", 1.21, 1.84), "Si": ("#F0C8A0", 1.20, 2.10),
    "P": ("#FF8000", 1.05, 1.80), "S": ("#FFFF30", 1.02, 1.80),
    "Cl": ("#1FF01F", 0.99, 1.75), "Ar": ("#80D1E3", 1.51, 1.88),
    "K": ("#8F40D4", 2.03, 2.75), "Ca": ("#3DFF00", 1.76, 2.31),
    "Sc": ("#E6E6E6", 1.70, 2.11), "Ti": ("#BFC2C7", 1.60, 2.00),
    "V": ("#A6A6AB", 1.53, 2.00), "Cr": ("#8A99C7", 1.39, 2.00),
    "Mn": ("#9C7AC7", 1.61, 2.00), "Fe": ("#E06633", 1.52, 2.00),
    "Co": ("#F090A0", 1.26, 2.00), "Ni": ("#50D050", 1.24, 1.63),
    "Cu": ("#C88033", 1.32, 1.40), "Zn": ("#7D80B0", 1.22, 1.39),
    "Ga": ("#C28F8F", 1.22, 1.87), "Ge": ("#668F8F", 1.17, 2.11),
    "As": ("#BD80E3", 1.21, 1.85), "Se": ("#FFA100", 1.22, 1.90),
    "Br": ("#A62929", 1.21, 1.85), "Kr": ("#5CB8D1", 1.50, 2.02),
    "Rb": ("#702EB0", 2.20, 3.03), "Sr": ("#00FF00", 1.95, 2.49),
    "Y": ("#94FFFF", 1.90, 2.00), "Zr": ("#94E0E0", 1.75, 2.00),
    "Nb": ("#73C2C9", 1.64, 2.00), "Mo": ("#54B5B5", 1.54, 2.00),
    "Tc": ("#3B9E9E", 1.47, 2.00), "Ru": ("#248F8F", 1.46, 2.00),
    "Rh": ("#0A7D8C", 1.42, 2.00), "Pd": ("#006985", 1.39, 1.63),
    "Ag": ("#C0C0C0", 1.45, 1.72), "Cd": ("#FFD98F", 1.44, 1.58),
    "In": ("#A67573", 1.42, 1.93), "Sn": ("#668080", 1.39, 2.17),
    "Sb": ("#9E63B5", 1.39, 2.06), "Te": ("#D47A00", 1.47, 2.06),
    "I": ("#940094", 1.40, 1.98), "Xe": ("#429EB0", 1.50, 2.16),
    "Cs": ("#57178F", 2.44, 3.43), "Ba": ("#00C900", 2.15, 2.68),
    "La": ("#70D4FF", 2.07, 2.00), "Ce": ("#FFFFC7", 2.04, 2.00),
    "Pr": ("#D9FFC7", 2.03, 2.00), "Nd": ("#C7FFC7", 2.01, 2.00),
    "Sm": ("#8FFFC7", 1.98, 2.00), "Eu": ("#61FFC7", 1.98, 2.00),
    "Gd": ("#45FFC7", 1.96, 2.00), "Tb": ("#30FFC7", 1.94, 2.00),
    "Dy": ("#1FFFC7", 1.92, 2.00), "Ho": ("#00FF9C", 1.92, 2.00),
    "Er": ("#00E675", 1.89, 2.00), "Tm": ("#00D452", 1.90, 2.00),
    "Yb": ("#00BF38", 1.87, 2.00), "Lu": ("#00AB24", 1.87, 2.00),
    "Hf": ("#4DC2FF", 1.75, 2.00), "Ta": ("#4DA6FF", 1.70, 2.00),
    "W": ("#2194D6", 1.62, 2.00), "Re": ("#267DAB", 1.51, 2.00),
    "Os": ("#266696", 1.44, 2.00), "Ir": ("#175487", 1.41, 2.00),
    "Pt": ("#D0D0E0", 1.36, 1.75), "Au": ("#FFD123", 1.36, 1.66),
    "Hg": ("#B8B8D0", 1.32, 1.55), "Tl": ("#A6544D", 1.45, 1.96),
    "Pb": ("#575961", 1.46, 2.02), "Bi": ("#9E4FB5", 1.48, 2.07),
    "U": ("#008FFF", 1.96, 1.86),
}

DEFAULT = ("#FF1493", 1.50, 2.00)

# metals that may act as a coordination centre (transition metals, lanthanides, actinides, some main-group metals)
METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Al", "Ga", "In",
    "Tl", "Sn", "Pb", "Bi",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "U",
}
# d-block metals, preferred as the centre in spin-crossover work
TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}


def is_element(sym: str) -> bool:
    return sym in _DATA


def color(sym: str) -> str:
    return _DATA.get(sym, DEFAULT)[0]


def covalent_radius(sym: str) -> float:
    return _DATA.get(sym, DEFAULT)[1]


def vdw_radius(sym: str) -> float:
    return _DATA.get(sym, DEFAULT)[2]


def normalize_symbol(text: str) -> str:
    """Extract the element symbol from strings such as 'Fe2+', 'FE', 'Fe1A' or 'N12'."""
    letters = ""
    for ch in text.strip():
        if ch.isalpha():
            letters += ch
        else:
            break
    if not letters:
        return "X"
    two = letters[:2].capitalize()
    one = letters[:1].upper()
    if len(letters) >= 2 and two in _DATA:
        # 'CL1' -> Cl; labels such as 'CA' (C-alpha) are rare in small-molecule CIFs, so treat as an element
        return two
    if one in _DATA:
        return one
    return two if two else "X"


# standard atomic weights (g/mol)
MASSES = {
    "H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Si": 28.085,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904, "Se": 78.971, "As": 74.922,
    "Li": 6.94, "Na": 22.990, "K": 39.098, "Mg": 24.305, "Ca": 40.078, "Al": 26.982,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845,
    "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630,
    "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07,
    "Rh": 102.906, "Pd": 106.42, "Ag": 107.868, "Cd": 112.414, "In": 114.818, "Sn": 118.710,
    "Sb": 121.760, "Te": 127.60, "La": 138.905, "Ce": 140.116, "Pr": 140.908, "Nd": 144.242,
    "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.925, "Dy": 162.500, "Ho": 164.930,
    "Er": 167.259, "Tm": 168.934, "Yb": 173.045, "Lu": 174.967, "Hf": 178.49, "Ta": 180.948,
    "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084, "Au": 196.967,
    "Hg": 200.592, "Tl": 204.38, "Pb": 207.2, "Bi": 208.980, "U": 238.029,
}


def mass(sym: str) -> float:
    return MASSES.get(sym, 50.0)
