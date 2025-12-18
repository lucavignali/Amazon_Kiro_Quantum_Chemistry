#!/usr/bin/env python3
"""
Generate PyMOL script to load any combination of cube files.
Automatically creates unique object names based on filenames.

Usage: python load_cubes.py H_orbital_1*.cube O_orbital_4*.cube H2O_orbital_3*.cube
"""
import sys
import glob
import os

def generate_pml(cube_files):
    """Generate PyMOL script to load cube files with unique names."""
    
    # Expand wildcards
    expanded = []
    for pattern in cube_files:
        matches = glob.glob(pattern)
        if matches:
            expanded.extend(matches)
        else:
            expanded.append(pattern)
    
    # Find corresponding PDB files and extract molecule names
    pdb_files = set()
    mol_names = set()
    for cube in expanded:
        base = os.path.basename(cube)
        mol_name = base.split('_orbital_')[0]
        mol_names.add(mol_name)
        pdb = f"{mol_name}_structure.pdb"
        if os.path.exists(pdb):
            pdb_files.add(pdb)
    
    # Generate output filename from molecule names
    output = '_'.join(sorted(mol_names)) + '_compare.pml'
    
    with open(output, 'w') as f:
        f.write("# Auto-generated PyMOL script\n")
        f.write(f"# Comparing: {', '.join(sorted(mol_names))}\n\n")
        f.write("delete all\n")
        f.write("bg_color white\n\n")
        
        # Load structures
        for pdb in sorted(pdb_files):
            mol_name = pdb.replace('_structure.pdb', '')
            f.write(f"# {mol_name} structure\n")
            f.write(f"load {pdb}, {mol_name}_atoms\n")
            f.write(f"show spheres, {mol_name}_atoms\n")
            f.write(f"set sphere_scale, 0.3, {mol_name}_atoms\n\n")
        
        # Load cubes
        colors = ['marine', 'forest', 'orange', 'purple', 'yellow', 'cyan', 'magenta', 'salmon']
        for i, cube in enumerate(sorted(expanded)):
            base = os.path.basename(cube).replace('.cube', '')
            color = colors[i % len(colors)]
            
            f.write(f"# {base}\n")
            f.write(f"load {cube}\n")
            f.write(f"isosurface {base}_pos, {base}, 0.05\n")
            f.write(f"isosurface {base}_neg, {base}, -0.05\n")
            f.write(f"color {color}, {base}_pos\n")
            f.write(f"color gray, {base}_neg\n")
            f.write(f"set transparency, 0.3, {base}_*\n\n")
        
        f.write("# Rendering\n")
        f.write("set ray_trace_mode, 1\n")
        f.write("set ray_shadows, 0\n")
        f.write("set antialias, 2\n")
        f.write("zoom all\n")
    
    print(f"Generated: {output}")
    print(f"Loaded {len(expanded)} cube files from {len(mol_names)} molecules: {', '.join(sorted(mol_names))}")
    print(f"\nTo visualize:")
    print(f"  pymol {output}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python load_cubes.py <cube_files...>")
        print("\nExamples:")
        print("  python load_cubes.py H_orbital_1*.cube O_orbital_4*.cube")
        print("  python load_cubes.py *_orbital_1_*.cube")
        print("  python load_cubes.py H2O_orbital_*.cube")
        sys.exit(1)
    
    generate_pml(sys.argv[1:])
