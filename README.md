# Virtual Screening Workflows: AutoDock Vina and DiffDock

This repository contains two complementary workflows demonstrating both **classical** and **AI-based** molecular docking approaches for structure-based virtual screening.  
The goal of this project is to establish an end-to-end reproducible pipeline — from protein and ligand preparation to post-docking analysis — using two paradigms:

- **AutoDock Vina:** Traditional physics-based scoring and docking  
- **DiffDock:** Deep learning–driven generative docking via diffusion models

---

## Repository Structure

```
autodock-workflow/      # Classical docking workflow using AutoDock Vina
DiffDock-workflow/      # AI-based docking workflow using DiffDock
```

---

## Workflow Overview

| Workflow | Tool | Description |
|-----------|------|-------------|
| **AutoDock Vina** | Physics-based | End-to-end molecular docking pipeline including receptor/ligand preparation, docking with AutoDock Vina, post-processing, ADME triage, and clustering analysis. |
| **DiffDock** | Machine learning–based | Equivalent workflow built using DiffDock for fast and accurate docking via diffusion models, including setup scripts, CSV generation, and post-docking clustering analysis. |

Both workflows use **CDK2 (PDB ID: 1H1Q)** as a model system, docking analogs of the co-crystal ligand *2A6 (NU2058)*.  
Ligands are sourced from a subset of the **Enamine Hinge Binder Library**.

---

## Environment Setup

Due to compatibility differences between the toolchains, the workflows are currently maintained in **separate conda environments**:

- `autodock_env` — includes AutoDock Vina, Meeko, Open Babel, and RDKit  
- `diffdock_env` — includes DiffDock, PyTorch, RDKit, Open Babel, and torch_geometric
- **Overall** :
     - DiffDock's poses were much reliable, however some cyclic ring torsions could be improved using GNINA minimization
     - AutoDock's rigid docking even with a maximum exhaustiveness of 32 did not result in good poses --> AutoDock Flex could be used or UniDock/GNINA/SMINA could be used instead
          - Other alternatives listed in the summary notebook and the attached summary slide.

```bash
# For AutoDock Vina
cd autodock-workflow
conda env create -f environment.yml
conda activate autodock_env

# For DiffDock
cd ../DiffDock-workflow
conda env create -f environment.yml
conda activate diffdock_env
```

**Note:**  
With additional dependency reconciliation (especially between PyTorch and RDKit versions), these workflows could ideally be merged into a **single unified environment** to streamline reproducibility.

---

## Suggested Execution Order

1. **AutoDock Workflow (`autodock-workflow/`)**
   - Prepare receptor and ligand libraries.
   - Perform docking using AutoDock Vina.
   - Analyze poses and cluster results.

2. **DiffDock Workflow (`DiffDock-workflow/`)**
   - Reuse same receptor and ligand datasets.
   - Prepare DiffDock input CSVs and run inference.
   - Analyze DiffDock poses and compare with AutoDock outputs.

---

## Outputs

- **AutoDock:** Top poses (SDF/CSV), ranked binding energies, clustering and ADME triage notebooks.  
- **DiffDock:** AI-predicted docking poses with confidence scores and ranked summaries.  

The combined analysis provides a side-by-side comparison of **physics-based vs AI-based docking performance** on the same target and ligand set.

---

## Summary
- This is included in the attached presentation **Summary-CADD-VS-Workflow**
- But also for the individual detailed summaries, two notebook files have been added, in the individual folders :
        - `step5-autodock_summary_selection.ipynb`
        - `4_diffdock-summary-analysis.ipynb`


## Citation

If referencing this work, please cite:

- Trott & Olson, *AutoDock Vina: Improving the speed and accuracy of docking*, *J. Comput. Chem.* (2010)  
- Corso et al., *DiffDock: Diffusion steps, twists, and turns for molecular docking*, *NeurIPS* (2022)  
- O’Boyle et al., *Open Babel: An open chemical toolbox*, *J. Cheminf.* (2011)
