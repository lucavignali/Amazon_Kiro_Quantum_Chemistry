# Quantum Orbital Calculator & Visualizer

A Python toolkit for computing and visualizing molecular orbitals using high-accuracy quantum chemistry methods. Built on PySCF, this project makes it easy to explore the electronic structure of atoms and molecules through 3D visualizations.

## What This Does

This toolkit bridges quantum chemistry calculations and visual exploration:

- **Calculate** molecular orbitals using state-of-the-art methods (CCSD - Coupled Cluster Singles Doubles)
- **Export** wavefunctions to standard cube file format for visualization
- **Visualize** orbitals in PyMOL with automatic isosurface rendering
- **Compare** orbitals across different molecules side-by-side

Perfect for educational purposes, research visualization, or understanding molecular electronic structure.

## Features

- Pre-configured molecules (atoms, diatomics, polyatomics) ready to calculate
- High-accuracy CCSD calculations with correlation-consistent basis sets
- Automatic cube file generation for any orbital
- PyMOL integration with customizable visualization scripts
- Multi-molecule comparison tools
- Orbital energy and occupancy analysis

## Prerequisites

**Software:**
- Python 3.8 or higher
- PyMOL (for visualization)

**Python Libraries:**
```bash
pip install -r requirements.txt
```
- PySCF >= 2.0.0 (quantum chemistry calculations)
- h5py >= 3.0.0 (data storage)
- numpy >= 1.20.0 (numerical operations)

**PyMOL Installation:**
- macOS: `brew install pymol`
- Linux: `sudo apt-get install pymol` or `conda install -c conda-forge pymol-open-source`
- Windows: Download from https://pymol.org/

See `INSTALL.md` for detailed setup instructions.

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
