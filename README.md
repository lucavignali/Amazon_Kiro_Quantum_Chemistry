# Quantum Orbital Calculator & Visualizer

Calculate and visualize molecular orbitals using PySCF with the most accurate methods (CCSD).

## Quick Start

### 1. Calculate Orbitals (Most Accurate Method)

```bash
# Atoms
python main.py H CCSD cc-pvtz
python main.py O CCSD cc-pvtz

# Molecules
python main.py H2O CCSD cc-pvtz
python main.py HCOOH CCSD cc-pvtz
```

**Output:** `results_<molecule>_CCSD_cc-pvtz_<timestamp>.h5`

**Method:** CCSD (Coupled Cluster Singles Doubles) - gold standard for molecular calculations
**Basis:** cc-pvtz (correlation-consistent polarized valence triple-zeta) - high accuracy

### 2. Export Orbitals to Cube Files

```bash
# Export all occupied orbitals
python export_to_pymol.py results_H2O_CCSD_cc-pvtz_*.h5

# Export specific orbitals (0-indexed)
python export_to_pymol.py results_O_CCSD_cc-pvtz_*.h5 --orbitals 0,3,4
```

**Output:**
- `<molecule>_structure.pdb` - atomic positions
- `<molecule>_orbital_N_E<energy>eV.cube` - orbital wavefunction (one per orbital)
- `<molecule>_pymol.pml` - PyMOL visualization script

**Note:** Each cube file contains one orbital as a 3D grid of wavefunction values.

### 3. Visualize in PyMOL

#### Single Molecule
```bash
pymol H2O_pymol.pml
```

#### Compare Multiple Molecules
```bash
# Compare H, O, and H2O orbitals
python load_cubes.py H_orbital_1*.cube O_orbital_4*.cube H2O_orbital_3*.cube

# Output: H_H2O_O_compare.pml
pymol H_H2O_O_compare.pml
```

**Features:**
- Auto-generates comparison scripts from any cube files
- Unique object names (no collisions)
- Supports wildcards: `*_orbital_1_*.cube`
- Output filename based on molecules included

## Available Molecules

Predefined in `molecules.py`:
- **Atoms:** H, He, Li, Be, B, C, N, O, F, Ne, Cl
- **Diatomics:** H2, Li2, B2, C2, N2, O2, F2, HF, LiH, OH, CO, NO
- **Polyatomics:** H2O, NH3, CH4, CO2, H2O2, HCl, H2SO4, HCOOH

## Key Files

- `main.py` - Orbital calculator
- `export_to_pymol.py` - Create cube files for visualization
- `load_cubes.py` - Generate PyMOL comparison scripts
- `molecules.py` - Molecular definitions

## Understanding Results

**Orbital energies:**
- More negative = more stable/bound
- Example: O atom 1s (-27.2 eV) vs 2p (-13.5 eV)

**Occupancy:**
- 1.0 = fully occupied (2 electrons)
- 0.5 = half-occupied (1 electron, unpaired)
- 0.0 = virtual/unoccupied orbital

**Isosurfaces in PyMOL:**
- Positive (colored): ψ > 0
- Negative (gray): ψ < 0
- Threshold typically 0.02-0.05

## Full Documentation

See `README_full.md` for complete documentation including:
- All calculation methods (HF, DFT, CCSD)
- Geometry optimization
- Superposition animations
- Unreal Engine export
- Detailed theory and examples
