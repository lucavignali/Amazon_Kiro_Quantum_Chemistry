import sys
import numpy as np
import scipy.linalg
from pyscf import gto, scf, dft, cc
from pyscf.geomopt.geometric_solver import optimize
from datetime import datetime
import h5py
from molecules import MOLECULES

def analyze_orbital_character(mol, mo_coeff, orbital_idx):
    """Analyze atomic contributions to a molecular orbital using Mulliken population."""
    # Handle different mo_coeff formats
    if isinstance(mo_coeff, tuple):
        # Unrestricted: (alpha, beta) tuple
        mo = mo_coeff[0][:, orbital_idx]
    elif mo_coeff.ndim == 3:
        # Unrestricted: (2, nao, nmo) array
        mo = mo_coeff[0][:, orbital_idx]
    elif mo_coeff.ndim == 2:
        # Restricted: (nao, nmo) array
        mo = mo_coeff[:, orbital_idx]
    elif mo_coeff.ndim == 1:
        # Single orbital: (nao,) array
        mo = mo_coeff
    else:
        # Fallback
        return "Unknown"
    
    ao_labels = mol.ao_labels(fmt=None)
    
    atom_contributions = {}
    for i in range(mol.natm):
        atom_symbol = mol.atom_symbol(i)
        atom_key = f"{atom_symbol}{i+1}"
        contrib = sum(mo[ao_idx]**2 for ao_idx, (atom_idx, _, _, _) in enumerate(ao_labels) if atom_idx == i)
        atom_contributions[atom_key] = contrib
    
    total = sum(atom_contributions.values())
    if total > 0:
        atom_contributions = {k: v/total*100 for k, v in atom_contributions.items()}
    
    max_contrib = max(atom_contributions.values()) if atom_contributions else 0
    if max_contrib > 85:
        return "Non-bonding"
    elif max_contrib < 65:
        return "Bonding"
    else:
        return "Weak-bond"

def check_stability(mf, max_cycles=5):
    """Check and optimize wavefunction stability"""
    for cycle in range(max_cycles):
        mo1, _, stable, _ = mf.stability(return_status=True)
        if stable:
            print(f"✓ Wavefunction is stable (cycle {cycle})")
            return mf
        print(f"⚠ Instability detected, re-optimizing (cycle {cycle})")
        dm1 = mf.make_rdm1(mo1, mf.mo_occ)
        mf = mf.run(dm1)
    print("⚠ Stability optimization did not converge")
    return mf

def solve_molecule(molecule_name, method='DFT', basis='6-31g', xc='B3LYP', optimize_geom=False):
    """
    Solve multi-electron Schrödinger equation for diatomic molecules.
    
    Args:
        molecule_name: 'H2', 'O2', 'N2', 'F2', 'H2+', 'Li2', etc.
        method: 'HF' (Hartree-Fock), 'DFT' (Density Functional Theory), or 'CCSD' (Coupled Cluster)
        basis: Basis set (e.g., '6-31g', 'cc-pvdz')
        xc: Exchange-correlation functional for DFT (e.g., 'B3LYP', 'PBE')
        optimize_geom: If True, optimize geometry before calculating orbitals
    """
    
    if molecule_name not in MOLECULES:
        print(f"Unknown molecule: {molecule_name}")
        print(f"Available: {list(MOLECULES.keys())}")
        sys.exit(1)
    
    atom_string, charge, spin = MOLECULES[molecule_name]
    
    # Build molecule
    mol = gto.M(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        unit='Angstrom'
    )
    
    print(f"\n{'='*60}")
    print(f"Molecule: {molecule_name}")
    print(f"Method: {method}")
    print(f"Basis: {basis}")
    if method == 'DFT':
        print(f"XC Functional: {xc}")
    print(f"Charge: {charge}")
    print(f"Spin: {mol.spin}")
    print(f"Number of electrons: {mol.nelectron}")
    print(f"Number of basis functions: {mol.nao}")
    if optimize_geom:
        print(f"Geometry optimization: ENABLED")
    print(f"{'='*60}\n")
    
    # Optimize geometry if requested
    if optimize_geom:
        print("Optimizing molecular geometry...")
        # Use DFT for optimization (faster than CCSD)
        if mol.spin == 0:
            mf_opt = dft.RKS(mol)
        else:
            mf_opt = dft.ROKS(mol)
        mf_opt.xc = 'B3LYP'
        
        mol_opt = optimize(mf_opt, maxsteps=50)
        
        print("\nOptimized geometry:")
        for i in range(mol_opt.natm):
            coord = mol_opt.atom_coord(i)
            print(f"  {mol_opt.atom_symbol(i)}: {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f} Angstrom")
        
        # Calculate bond lengths
        if mol_opt.natm == 2:
            bond_length = np.linalg.norm(mol_opt.atom_coord(0) - mol_opt.atom_coord(1))
            print(f"\nOptimized bond length: {bond_length:.6f} Angstrom")
        
        mol = mol_opt
        print()
    
    # Solve
    if method == 'HF':
        if mol.spin == 0:
            mf = scf.RHF(mol)
        else:
            mf = scf.ROHF(mol)
    elif method == 'CCSD':
        # CCSD requires HF reference
        if mol.spin == 0:
            mf = scf.RHF(mol)
        else:
            mf = scf.ROHF(mol)
    else:  # DFT
        if mol.spin == 0:
            mf = dft.RKS(mol)
        else:
            mf = dft.ROKS(mol)
        mf.xc = xc
    
    energy = mf.kernel()
    
    # Check wavefunction stability
    print("\nChecking wavefunction stability...")
    mf = check_stability(mf)
    energy = mf.e_tot
    
    # Run CCSD if requested
    if method == 'CCSD':
        print("\nRunning CCSD calculation...")
        if mol.spin == 0:
            mycc = cc.CCSD(mf)
        else:
            mycc = cc.UCCSD(mf)
        mycc.kernel()
        energy = mycc.e_tot
        
        # Get natural orbitals from CCSD density matrix
        print("Computing natural orbitals...")
        dm1 = mycc.make_rdm1()  # 1-electron density matrix in MO basis
        
        if mol.spin == 0:
            # Diagonalize density matrix to get natural orbitals
            no_occ, no_coeff_mo = np.linalg.eigh(dm1)
            # Sort by occupancy (descending)
            idx = np.argsort(no_occ)[::-1]
            mo_occ = no_occ[idx]
            no_coeff_mo = no_coeff_mo[:, idx]
            # Transform to AO basis: C_AO = C_HF @ C_NO
            mo_coeff = mf.mo_coeff @ no_coeff_mo
            # Use negative occupancy as proxy for energy ordering
            mo_energy = -mo_occ
        else:
            # For unrestricted CCSD (UCCSD uses ROHF reference)
            # ROHF mo_coeff is already 2D: (nao, nmo), not a tuple
            dm1_avg = (dm1[0] + dm1[1]) / 2
            no_occ, no_coeff_mo = np.linalg.eigh(dm1_avg)
            idx = np.argsort(no_occ)[::-1]
            mo_occ = no_occ[idx]
            no_coeff_mo = no_coeff_mo[:, idx]
            # mf.mo_coeff is already (nao, nmo) for ROHF
            mo_coeff = mf.mo_coeff @ no_coeff_mo
            mo_energy = -mo_occ
    else:
        mo_energy = mf.mo_energy
        mo_occ = mf.mo_occ
        mo_coeff = mf.mo_coeff
        
        # Debug: check what we have
        print(f"DEBUG: mo_coeff type: {type(mo_coeff)}, shape: {mo_coeff.shape if hasattr(mo_coeff, 'shape') else 'N/A'}")
        
        # For unrestricted (open-shell) HF/DFT, extract alpha orbitals
        # ROHF returns 2D array directly, UHF/UOKS returns tuple
        if isinstance(mo_coeff, tuple):
            mo_energy = mo_energy[0]
            mo_occ = mo_occ[0]
            mo_coeff = mo_coeff[0]  # Extract alpha: (nao, nmo)
        # mo_coeff should already be 2D for ROHF
    
    # Display results
    print(f"\nTotal Energy: {energy:.6f} Hartree ({energy * 27.211:.2f} eV)")
    print(f"\nMolecular Orbital Energies (eV):")
    print(f"{'Orbital':<10} {'Energy (eV)':<15} {'Occupancy':<12} {'Character':<15}")
    print("-" * 55)
    
    for i, (e, occ) in enumerate(zip(mo_energy, mo_occ)):
        orb_char = analyze_orbital_character(mol, mo_coeff, i)
        print(f"{i+1:<10} {e*27.211:<15.2f} {occ:<12.1f} {orb_char:<15}")
    
    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    opt_suffix = "_opt" if optimize_geom else ""
    filename = f"results_{molecule_name}_{method}_{basis}{opt_suffix}_{timestamp}.h5"
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('mo_energy', data=mo_energy)
        f.create_dataset('mo_occ', data=mo_occ)
        f.create_dataset('mo_coeff', data=mo_coeff)
        f.create_dataset('total_energy', data=energy)
        f.attrs['molecule'] = molecule_name
        f.attrs['method'] = method
        f.attrs['basis'] = basis
        f.attrs['nelectron'] = mol.nelectron
        f.attrs['spin'] = mol.spin
        
        # Save molecular geometry
        coords = mol.atom_coords()
        f.create_dataset('atom_coords', data=coords)
        f.attrs['atom_symbols'] = [mol.atom_symbol(i) for i in range(mol.natm)]
    
    print(f"\nResults saved to {filename}")
    return filename

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py <molecule> [method] [basis] [xc] [--optimize]")
        print("\nExamples:")
        print("  python main.py H2")
        print("  python main.py O2 DFT cc-pvdz")
        print("  python main.py H2+ HF")
        print("  python main.py H2O CCSD cc-pvtz --optimize")
        sys.exit(1)
    
    # Check for --optimize flag
    optimize_geom = '--optimize' in sys.argv
    args = [arg for arg in sys.argv[1:] if arg != '--optimize']
    
    molecule = args[0]
    method = args[1] if len(args) > 1 else 'DFT'
    basis = args[2] if len(args) > 2 else '6-31g'
    xc = args[3] if len(args) > 3 else 'B3LYP'
    
    solve_molecule(molecule, method, basis, xc, optimize_geom)
