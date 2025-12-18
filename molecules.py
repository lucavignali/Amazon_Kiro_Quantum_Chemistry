"""
Molecular definitions for diatomic molecules.
Each entry: (atom_string, charge, spin)
"""

MOLECULES = {
    'H2': ('H 0 0 0; H 0 0 0.74', 0, 0),
    'H2+': ('H 0 0 0; H 0 0 1.06', 1, 0),
    'H': ('H 0 0 0', 0, 1),  # Hydrogen atom
    'O': ('O 0 0 0', 0, 2),  # Oxygen atom (triplet ground state)
    'Cl': ('Cl 0 0 0', 0, 1),  # Chlorine atom (doublet ground state)
    'Li2': ('Li 0 0 0; Li 0 0 2.67', 0, 0),
    'B2': ('B 0 0 0; B 0 0 1.59', 0, 2),
    'C2': ('C 0 0 0; C 0 0 1.24', 0, 0),
    'N2': ('N 0 0 0; N 0 0 1.10', 0, 0),
    'O2': ('O 0 0 0; O 0 0 1.21', 0, 2),
    'F2': ('F 0 0 0; F 0 0 1.41', 0, 0),
    'HF': ('H 0 0 0; F 0 0 0.92', 0, 0),
    'HCl': ('H 0 0 0; Cl 0 0 1.27', 0, 0),
    'CO': ('C 0 0 0; O 0 0 1.13', 0, 0),
    'H2O': ('O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0', 0, 0),
    'HCOOH': ('C 0 0 0; O 1.202 0 0; O -0.544 1.188 0; H -0.544 -1.102 0; H -1.489 1.188 0', 0, 0),
    'H2SO4': ('S 0 0 0; O 1.422 0 0; O -0.711 1.232 0; O -0.711 -1.232 0; O 0 0 1.574; H -1.661 1.232 0; H -1.661 -1.232 0', 0, 0),
    'LiCoO2': ('Co 0 0 0; O 0 0 1.93; O 0 0 -1.93; O 1.93 0 0; O -1.93 0 0; O 0 1.93 0; O 0 -1.93 0; Li 0 0 4.7', 0, 0),  # CoO6 octahedron + Li
    'H_O_far': ('H 0 0 0; O 0 0 10.0', 0, 3),  # Non-interacting H and O atoms (10 Å separation)
}
