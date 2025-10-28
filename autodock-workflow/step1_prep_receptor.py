#!/usr/bin/env python3
"""
step1_prep_receptor.py
------------------------
Prepares the CDK2 (1H1Q) receptor for docking.

Steps:
1. Download PDB 1H1Q
2. Select protein chains A & B (CDK2–Cyclin A complex)
3. Remove water and heteroatoms
4. Determine docking box center from co-ligand (2A6)
5. Prepare receptor with Meeko (PDBQT + box)

Usage:
    python step1_prep_receptor.py --pdb_id 1H1Q --ligand_resname 2A6 --chains A B

Requirements:
    conda install -c conda-forge prody meeko
"""

import os
import time
import argparse
import subprocess
from prody import parsePDB, writePDB, calcCenter

def main():
    parser = argparse.ArgumentParser(
        description="Prepare CDK2 (1H1Q) receptor for AutoDock Vina or DiffDock"
    )
    parser.add_argument("--pdb_id", default="1H1Q", help="RCSB PDB ID")
    parser.add_argument("--ligand_resname", default="2A6", help="Ligand residue name (e.g., 2A6)")
    parser.add_argument("--box_size", nargs=3, type=float, default=[20.0, 20.0, 20.0],
                        help="Docking box dimensions in Å (x y z)")
    parser.add_argument("--chains", nargs="+", default=["A", "B"], help="Protein chain IDs to keep")
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    ligand_resname = args.ligand_resname
    chains = args.chains
    box_size = args.box_size

    print(f"[INFO] Downloading {pdb_id} from RCSB...")
    os.system(f"curl -s https://files.rcsb.org/view/{pdb_id}.pdb -o {pdb_id}.pdb")

    # --- Step 1: Parse PDB ---
    print(f"[INFO] Parsing {pdb_id}.pdb...")
    atoms = parsePDB(pdb_id)
    if atoms is None:
        raise FileNotFoundError(f"Failed to parse {pdb_id}.pdb")

    # --- Step 2: Select chains A & B, remove water & heteroatoms ---
    chain_selection = " or ".join([f"chain {c}" for c in chains])
    print(f"[INFO] Selecting receptor atoms ({', '.join(chains)}), excluding water/heteroatoms...")
    receptor_atoms = atoms.select(f"({chain_selection}) and not water and not hetero")
    if receptor_atoms is None:
        raise ValueError("No receptor atoms found. Check chain IDs or PDB content.")

    pdb_path = f"{pdb_id}_receptor_atoms.pdb"
    writePDB(pdb_path, receptor_atoms)

    # --- Wait for file creation ---
    for _ in range(10):
        if os.path.exists(pdb_path):
            break
        time.sleep(0.2)
    else:
        raise FileNotFoundError(f"{pdb_path} not found after writing.")

    # --- Step 3: Calculate docking box center ---
    print(f"[INFO] Calculating docking box center using co-ligand: {ligand_resname}")
    ligand_atoms = atoms.select(f"chain {chains[0]} and resname {ligand_resname}")
    if ligand_atoms is None:
        raise ValueError(f"Ligand '{ligand_resname}' not found in chain {chains[0]}.")

    cx, cy, cz = calcCenter(ligand_atoms)
    print(f"[INFO] Box center = ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    # --- Step 4: Prepare receptor with Meeko ---
    renamed_pdb = f"receptor_{pdb_id}.pdb"
    os.rename(pdb_path, renamed_pdb)
    out_base = f"{pdb_id}_receptorH"

    print("[INFO] Running Meeko receptor preparation...")
    cmd = [
        "mk_prepare_receptor.py",
        "-i", renamed_pdb,
        "-o", out_base,
        "-p", "-v",
        "--box_center", str(cx), str(cy), str(cz),
        "--box_size", str(box_size[0]), str(box_size[1]), str(box_size[2])
    ]
    subprocess.run(cmd, check=True)
    print(f"[SUCCESS] Receptor prepared → {out_base}.pdbqt and {out_base}.box.pdb")

    # --- Step 5: Print summary ---
    print("\n=== Summary ===")
    print(f"PDB ID:          {pdb_id}")
    print(f"Chains retained: {', '.join(chains)}")
    print(f"Ligand:          {ligand_resname}")
    print(f"Box center:      ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"Box size:        {box_size}")
    print(f"Output files:    {out_base}.pdbqt, {out_base}.box.pdb")
    print("================")

if __name__ == "__main__":
    main()

