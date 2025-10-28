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

- Both workflows use **CDK2 (PDB ID: 1H1Q)** as a model system, docking analogs of the co-crystal ligand *2A6 (NU2058)*.  
- Ligands are sourced from a subset of the **Enamine Hinge Binder Library** and close-in analogs (80% similarity) with co-crystal ligand 2A6, from PubChem.
   - the co-crystal ligand was also used in the workflow, for a quic validation of the process - although not include in hit triaging/candate selection.

---

## Environment Setup

Due to compatibility differences between the toolchains, the workflows are currently maintained in **separate conda environments**:

- `autodock_env` — includes AutoDock Vina, Meeko, Open Babel, and RDKit  
- `diffdock_env` — includes DiffDock, PyTorch, RDKit, Open Babel, and torch_geometric

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
   - Reuse same receptor and ligand datasets, by converting the formats (bash scripts included)
   - Prepare DiffDock input CSVs and run inference.
   - Analyze DiffDock poses and compare with AutoDock outputs.

---

## Outputs

- **AutoDock:** Top poses (SDF/CSV), ranked binding energies, clustering and ADME triage notebooks.  
- **DiffDock:** AI-predicted docking poses with confidence scores and ranked summaries.  

The combined analysis provides a side-by-side comparison of **physics-based vs AI-based docking performance** on the same target and ligand set.

---
## Final List of Candidates 
Autodock was completed and produced selections; however, due to low confidence in hinge-binder poses, I did not advance candidates from those. The list below reflects compounds selected via **DiffDock**.

### Selected Compounds

| Catalog_ID   | Source  | Method   |
|--------------|---------|----------|
| Z16226232    | Enamine | DiffDock |
| Z3071631534  | Enamine | DiffDock |
| 49863196     | PubChem | DiffDock |
| Z56854611    | Enamine | DiffDock |
| Z1013695804  | Enamine | DiffDock |
| Z1514091612  | Enamine | DiffDock |
| 5329551      | PubChem | DiffDock |
| Z2636563025  | Enamine | DiffDock |
| Z442313918   | Enamine | DiffDock |
| 5329506      | PubChem | DiffDock |

## Notes
- Autodock runs completed successfully, but hinge-binder pose confidence was insufficient to move those candidates forward.
- DiffDock selections prioritized pose plausibility near the hinge region.
- The PubChem compounds are already CDK2 inhibitors and part of patents, but the data could be explored for further evaluation, and possibly analysis on lead optimization.
  
## Citation

If referencing this work, please cite:

- Trott & Olson, *AutoDock Vina: Improving the speed and accuracy of docking*, *J. Comput. Chem.* (2010)  
- Corso et al., *DiffDock: Diffusion steps, twists, and turns for molecular docking*, *NeurIPS* (2022)  
- O’Boyle et al., *Open Babel: An open chemical toolbox*, *J. Cheminf.* (2011)
