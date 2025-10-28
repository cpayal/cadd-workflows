# AutoDock Vina Virtual Screening Workflow

This repository contains a minimal reproducible module for structure-based virtual screening using **AutoDock Vina**.  
It demonstrates ligand preparation, docking, and analysis using an example target (**CDK2**, PDB ID: *1H1Q*).

---

## Repository Structure

```
step1_prep_receptor.py            # Protein preprocessing: clean PDB, extract chains, add hydrogens
step2-prep-ligs.py                # Ligand preparation from SDF (sample workflow)
step3-batch_dock_vina.py          # Batch docking with AutoDock Vina using multiprocessing
example-step3-run.txt             # Example CLI command for Step 3
step4_process_docking_outputs.py  # Parse docking results, extract top poses, generate summary CSV/SDF
step5-autodock_summary_selection.ipynb  # Analysis and visualization of docking results
Enamine_Hinge_Binding_sample.sdf  # Example ligand library (subset of Enamine hinge binder set)
environment.yml                   # Conda environment file for reproducibility
example_outputs/                  # Example docking outputs (poses, logs, summaries)
notebooks/                        # Additional Jupyter notebooks for analysis
logs/                             # Run-time logs
test/                             # Test dataset and debug scripts
```

---

## Environment Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate vina_env
```

**Environment includes:**
- Python ≥ 3.9  
- RDKit for cheminformatics  
- Meeko and AutoDockTools for receptor/ligand preparation  
- Open Babel for format conversions  
- pandas, tqdm, numpy for data handling  

---

## Workflow Overview

| Step | Script | Description |
|------|---------|-------------|
| 1 | step1_prep_receptor.py | Prepares the protein for docking: removes water, isolates target chains, adds hydrogens and charges, saves PDBQT. |
| 2 | step2-prep-ligs.py | Converts ligands from SDF to PDBQT using RDKit + Meeko; filters duplicates. |
| 3 | step3-batch_dock_vina.py | Runs AutoDock Vina in parallel for all ligands; outputs best poses and energies. |
| 4 | step4_process_docking_outputs.py | Aggregates docking logs, extracts poses, calculates ΔG (kcal/mol), exports SDF and CSV summaries. |
| 5 | step5-autodock_summary_selection.ipynb | Jupyter notebook for hit selection, clustering (ECFP + Tanimoto + Butina), and visualization. |

---

## Example Run

```bash
# Step 1: Prepare receptor
python step1_prep_receptor.py --pdb 1H1Q.pdb --chains A,B --out receptor_prepared.pdbqt

# Step 2: Prepare ligands
python step2-prep-ligs.py --sdf Enamine_Hinge_Binding_sample.sdf --out ligands_pdbqt/

# Step 3: Run docking
python step3-batch_dock_vina.py --receptor receptor_prepared.pdbqt     --ligands 'ligands_pdbqt/*.pdbqt'     --outdir docked_out     --center_x 10.5 --center_y 15.2 --center_z 22.3     --size_x 20 --size_y 20 --size_z 20

# Step 4: Process results
python step4_process_docking_outputs.py --glob 'docked_out/*.pdbqt' --outdir exports
```

Run `step5-autodock_summary_selection.ipynb` to reproduce clustering and hit selection.

---

## Output Summary

- **Pose energies CSV:** Binding energies per ligand  
- **Hits SDF:** Top poses filtered by RMSD and ΔG  
- **Cluster plots:** ECFP-based similarity clusters visualized in the notebook  
- **Insights:** PubChem analogs clustered with the co-crystal ligand; diverse scaffolds selected for follow-up  

---

## Notes

- Example target: **CDK2 (PDB ID: 1H1Q)** with ligand **2A6 (NU2058)**  
- Ligand library: **Enamine Hinge Binder sample subset**  
- Top poses show ≤ 1 Å RMSD agreement between Vina and DiffDock predictions  

---

## Citation

If reusing this workflow, please cite:  
- Trott & Olson, *AutoDock Vina: Improving the speed and accuracy of docking*, *J. Comput. Chem.* (2010)  
- O’Boyle et al., *Open Babel: An open chemical toolbox*, *J. Cheminf.* (2011)
