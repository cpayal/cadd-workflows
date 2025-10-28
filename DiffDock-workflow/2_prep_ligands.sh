#!/usr/bin/env bash
# Usage: ./2_prep_ligands.sh <pdbqt_dir> <protein_pdb> <out_dir>
# Example: ./2_prep_ligands.sh ligands_pdbqt receptor.pdb results_dd
#
# Converts AutoDock PDBQT ligands to SDF and writes DiffDock CSV:
# complex_name,protein_path,ligand_description,protein_sequence,pocket_center_x,pocket_center_y,pocket_center_z,pocket_radius

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <pdbqt_dir> <protein_pdb> <out_dir>" >&2
  exit 1
fi

pdbqt_dir=$(realpath "$1")
protein_path=$(realpath "$2")
out_dir=$(realpath "$3")

# Pocket center and radius (as requested)
POCKET_X=6.234
POCKET_Y=44.214
POCKET_Z=50.819
POCKET_RADIUS=10.0

if [[ ! -d "$pdbqt_dir" ]]; then
  echo "[ERROR] PDBQT directory not found: $pdbqt_dir" >&2
  exit 1
fi
if [[ ! -f "$protein_path" ]]; then
  echo "[ERROR] Protein PDB not found: $protein_path" >&2
  exit 1
fi
if ! command -v obabel >/dev/null 2>&1; then
  echo "[ERROR] OpenBabel (obabel) not found. Try: conda install -c conda-forge openbabel" >&2
  exit 1
fi

mkdir -p "$out_dir/ligands_sdf"
csv="$out_dir/screen.csv"

echo "[INFO] Converting PDBQT → SDF..."
shopt -s nullglob
count=0
for f in "$pdbqt_dir"/*.pdbqt; do
  stem=$(basename "$f" .pdbqt)
  out="$out_dir/ligands_sdf/${stem}.sdf"
  # Add hydrogens; remove -h if you prefer original protonation
  obabel -ipdbqt "$f" -osdf -O "$out" -h >/dev/null 2>&1
  ((count++)) || true
done
if [[ $count -eq 0 ]]; then
  echo "[ERROR] No .pdbqt files found in $pdbqt_dir" >&2
  exit 1
fi

echo "[INFO] Generating CSV for DiffDock..."
# Header extended with pocket columns
echo "complex_name,protein_path,ligand_description,protein_sequence,pocket_center_x,pocket_center_y,pocket_center_z,pocket_radius" > "$csv"

for sdf in "$out_dir"/ligands_sdf/*.sdf; do
  name=$(basename "$sdf" .sdf)
  abs_sdf=$(realpath "$sdf")
  # Leave protein_sequence empty (… , , …) to avoid folding
  printf '%s,"%s","%s",,%s,%s,%s,%s\n' \
    "$name" "$protein_path" "$abs_sdf" \
    "$POCKET_X" "$POCKET_Y" "$POCKET_Z" "$POCKET_RADIUS" >> "$csv"
done

echo "[DONE] CSV ready: $csv"
