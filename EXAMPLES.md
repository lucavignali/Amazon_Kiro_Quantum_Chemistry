# Examples

## Basic Usage

### Single Atoms

```bash
# Hydrogen atom
python main.py H CCSD cc-pvtz
python export_to_pymol.py results_H_CCSD_cc-pvtz_*.h5
pymol H_pymol.pml

# Oxygen atom
python main.py O CCSD cc-pvtz
python export_to_pymol.py results_O_CCSD_cc-pvtz_*.h5
pymol O_pymol.pml
```

### Simple Molecules

```bash
# Water
python main.py H2O CCSD cc-pvtz
python export_to_pymol.py results_H2O_CCSD_cc-pvtz_*.h5
pymol H2O_pymol.pml

# Formic acid
python main.py HCOOH CCSD cc-pvtz
python export_to_pymol.py results_HCOOH_CCSD_cc-pvtz_*.h5
pymol HCOOH_pymol.pml
```

## Comparing Molecules

Compare hydrogen, oxygen, and water orbitals:

```bash
# Calculate all three
python main.py H CCSD cc-pvtz
python main.py O CCSD cc-pvtz
python main.py H2O CCSD cc-pvtz

# Export orbitals
python export_to_pymol.py results_H_CCSD_cc-pvtz_*.h5
python export_to_pymol.py results_O_CCSD_cc-pvtz_*.h5
python export_to_pymol.py results_H2O_CCSD_cc-pvtz_*.h5

# Create comparison visualization
python load_cubes.py H_orbital_1*.cube O_orbital_4*.cube H2O_orbital_3*.cube

# View comparison
pymol H_H2O_O_compare.pml
```

## Selecting Specific Orbitals

Export only specific orbitals (0-indexed):

```bash
# Export only HOMO and LUMO of water
python export_to_pymol.py results_H2O_CCSD_cc-pvtz_*.h5 --orbitals 4,5
```

## Battery Materials

```bash
# Lithium cobalt oxide (battery cathode)
python main.py LiCoO2 CCSD cc-pvtz
python export_to_pymol.py results_LiCoO2_CCSD_cc-pvtz_*.h5
pymol LiCoO2_pymol.pml
```

## Available Molecules

See `molecules.py` for the complete list:
- **Atoms:** H, O, Cl
- **Diatomics:** H2, N2, O2, CO, HF, HCl
- **Polyatomics:** H2O, NH3, CH4, CO2, HCOOH, H2SO4, LiCoO2
