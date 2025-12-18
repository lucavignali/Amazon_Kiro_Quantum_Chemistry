# Installation Guide

## Prerequisites

- Python 3.8 or higher
- PyMOL (for visualization)

## Setup

1. **Clone the repository:**
```bash
git clone <your-repo-url>
cd quantum-orbital-visualizer
```

2. **Create virtual environment:**
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. **Install dependencies:**
```bash
pip install -r requirements.txt
```

4. **Install PyMOL:**
- **macOS:** `brew install pymol`
- **Linux:** `sudo apt-get install pymol` or `conda install -c conda-forge pymol-open-source`
- **Windows:** Download from https://pymol.org/

## Quick Test

```bash
# Calculate water molecule orbitals
python main.py H2O CCSD cc-pvtz

# Export for visualization
python export_to_pymol.py results_H2O_CCSD_cc-pvtz_*.h5

# View in PyMOL
pymol H2O_pymol.pml
```

If you see the water molecule with orbital isosurfaces, you're all set!

## Troubleshooting

**PySCF installation fails:**
- Try: `pip install pyscf --no-cache-dir`
- Or use conda: `conda install -c pyscf pyscf`

**PyMOL not found:**
- Make sure PyMOL is in your PATH
- Try running `pymol` directly to verify installation

**Calculation takes too long:**
- Start with smaller molecules (H, H2)
- Use smaller basis sets (6-31g instead of cc-pvtz)
