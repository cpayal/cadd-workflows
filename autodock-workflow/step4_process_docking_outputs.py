#!/usr/bin/env python3
"""


step4 : Process AutoDock/Vina outputs with Meeko; emit only "hits" and a top-pose summary.

Definitions (per request)
-------------------------
- Convergence: complete-linkage cluster where ALL pairwise pose RMSDs < rmsd_thresh (Å).
- Hit: CID with at least 'hit_min_poses' poses in the converged cluster (default: 5).
- Viable poses: the poses that belong to the converged cluster for a hit.

Outputs
-------
- hits.sdf / hits.csv
    * Only CIDs that are hits
    * Only viable poses (i.e., poses in the converged cluster)
- top_pose_hits.sdf / top_pose_hits.csv
    * One row/record per hit (the best-scoring pose within the cluster)
    * Fields include SMILES and RMSD dispersion metrics:
        - rmsd_pairwise_std: stdev of all pairwise RMSDs within cluster
        - rmsd_to_best_std: stdev of RMSD to the best pose within cluster

Dependencies
------------
pip install meeko rdkit-pypi pandas numpy
(or conda-forge for rdkit; meeko via pip)
"""

import os
import sys
import glob
import math
import argparse
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdMolAlign, Descriptors, rdMolDescriptors, Crippen, Lipinski

# ---- Meeko import (fail fast with a helpful message) ----
try:
    from meeko import PDBQTMolecule, RDKitMolCreate
except Exception as e:
    sys.exit(
        "ERROR: Meeko is required (pip install meeko).\n"
        "This script relies on Meeko to rebuild RDKit molecules from Vina/AD-GPU outputs.\n"
        f"Details: {e}"
    )


# ============================ Utilities ============================

def stem_to_cid(stem: str) -> str:
    """Trim trailing '_dock' to get canonical CID from filename stem."""
    return stem[:-5] if stem.endswith("_dock") else stem


def parse_vina_scores_from_pdbqt(pdbqt_path: str) -> List[float]:
    """Parse 'REMARK VINA RESULT:' energies (Vina PDBQT)."""
    scores = []
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        scores.append(float(parts[3]))
                    except ValueError:
                        scores.append(float("nan"))
    return scores


def rdkit_physchem_props(mol: Chem.Mol) -> Dict[str, float]:
    """MWT, TPSA, logP, HBD, HBA; return NaNs if anything fails."""
    try:
        return {
            "MWT": Descriptors.MolWt(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "logP": Crippen.MolLogP(mol),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
        }
    except Exception:
        return {"MWT": math.nan, "TPSA": math.nan, "logP": math.nan, "HBD": math.nan, "HBA": math.nan}


def mol_to_smiles(mol: Chem.Mol) -> str:
    """Canonical isomeric SMILES from RDKit mol (connectivity only)."""
    try:
        m = Chem.Mol(mol)  # copy; coords don't matter for SMILES
        return Chem.MolToSmiles(m, isomericSmiles=True)
    except Exception:
        return ""


def trim_conformers(m: Chem.Mol, keep_n: int) -> Chem.Mol:
    """Keep only first 'keep_n' conformers (poses)."""
    n = m.GetNumConformers()
    if keep_n >= n:
        return m
    out = Chem.Mol(m)
    conf_ids = [c.GetId() for c in m.GetConformers()][:keep_n]
    out.RemoveAllConformers()
    for cid in conf_ids:
        out.AddConformer(Chem.Conformer(m.GetConformer(cid)), assignId=True)
    return out


def clone_with_conformer(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    """Return a 1-conformer copy; handy for SDF writing with per-pose props."""
    new_mol = Chem.Mol(mol)
    conf = mol.GetConformer(conf_id)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf), assignId=True)
    return new_mol


def rmsd_matrix(mol_with_confs: Chem.Mol) -> np.ndarray:
    """Pairwise RMSD (Å) across conformers; uses prbId/refId signature."""
    n = mol_with_confs.GetNumConformers()
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            val = rdMolAlign.GetBestRMS(mol_with_confs, mol_with_confs, prbId=j, refId=i)
            M[i, j] = M[j, i] = float(val)
    return M


def greedy_complete_linkage_cluster(rms: np.ndarray, scores: List[float], threshold: float) -> List[int]:
    """
    Largest set of pose indices where *every pair* < threshold Å (complete-linkage).
    Greedy: seed with best score, then add by score if within threshold to ALL current members.
    Ties on size → pick better mean score.
    """
    n = len(scores)
    order = sorted(range(n), key=lambda i: (scores[i] if not math.isnan(scores[i]) else float("inf")))
    best_cluster, best_mean = [], float("inf")
    for seed in order:
        cluster = [seed]
        for j in order:
            if j == seed:
                continue
            if all(rms[i, j] < threshold for i in cluster):
                cluster.append(j)
        if len(cluster) > len(best_cluster):
            best_cluster = cluster
            vals = [scores[i] for i in cluster if not math.isnan(scores[i])]
            best_mean = np.mean(vals) if vals else float("inf")
        elif len(cluster) == len(best_cluster) and cluster:
            vals = [scores[i] for i in cluster if not math.isnan(scores[i])]
            m = np.mean(vals) if vals else float("inf")
            if m < best_mean:
                best_cluster, best_mean = cluster, m
    return best_cluster


def read_pdbqt_with_meeko(path: str, max_poses: int = 10, is_dlg: bool = False) -> Tuple[Chem.Mol, List[float], List[int]]:
    """Meeko → RDKit (poses as conformers), plus Vina scores for PDBQT."""
    pdbqt_mol = PDBQTMolecule.from_file(path, is_dlg=is_dlg, skip_typing=True)
    rdkit_list = RDKitMolCreate.from_pdbqt_mol(
        pdbqt_mol,
        only_cluster_leads=False,
        keep_flexres=False
    )
    rdkit_list = [m for m in rdkit_list if m is not None]
    if not rdkit_list:
        raise ValueError("Meeko could not create RDKit molecules from this file.")
    mol = rdkit_list[0]
    if mol.GetNumConformers() > max_poses:
        mol = trim_conformers(mol, max_poses)

    if not is_dlg:
        scores = parse_vina_scores_from_pdbqt(path)
    else:
        scores = []  # add a DLG parser if you need it

    nconf = mol.GetNumConformers()
    if len(scores) < nconf:
        scores = scores + [float("nan")] * (nconf - len(scores))
    elif len(scores) > nconf:
        scores = scores[:nconf]

    pose_idxs = list(range(1, nconf + 1))
    return mol, scores, pose_idxs


def cluster_rmsd_stats(rms: np.ndarray, cluster_idx: List[int], best_idx: int) -> Dict[str, float]:
    """
    Compute dispersion metrics for a cluster:
      - rmsd_pairwise_std: stdev of all pairwise RMSDs within cluster (upper-tri)
      - rmsd_to_best_std: stdev of RMSDs from cluster poses to best pose (excl. self)
    """
    k = len(cluster_idx)
    if k <= 1:
        return {"rmsd_pairwise_std": float("nan"), "rmsd_to_best_std": float("nan")}

    sub = rms[np.ix_(cluster_idx, cluster_idx)]
    iu = np.triu_indices(k, k=1)
    pair_vals = sub[iu]
    rmsd_pairwise_std = float(np.std(pair_vals, ddof=1)) if pair_vals.size >= 2 else 0.0

    # best_idx is in original indexing; convert to cluster-local column
    best_local = cluster_idx.index(best_idx)
    to_best = np.delete(sub[best_local, :], best_local)  # drop self (0.0)
    rmsd_to_best_std = float(np.std(to_best, ddof=1)) if to_best.size >= 2 else (0.0 if to_best.size == 1 else float("nan"))

    return {"rmsd_pairwise_std": rmsd_pairwise_std, "rmsd_to_best_std": rmsd_to_best_std}


# ============================ Main ============================

def main():
    ap = argparse.ArgumentParser(description="Write only hits (≥ min poses in converged cluster) and top-pose summary.")
    ap.add_argument("--input_dir", required=True, help="Directory with .pdbqt (Vina) or .dlg (AD-GPU) files")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--rmsd_thresh", type=float, default=2.0, help="Convergence threshold Å (default 2.0)")
    ap.add_argument("--hit_min_poses", type=int, default=5, help="Minimum poses in converged cluster to call a hit (default 5)")
    ap.add_argument("--max_poses", type=int, default=10, help="Max poses per CID to consider (default 10)")
    ap.add_argument("--write-rmsd-matrices", action="store_true", help="Write per-CID pairwise RMSD CSVs")
    ap.add_argument("--smiles_csv", default="", help="Optional CSV (CID,smiles) to override SMILES calculation")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Gather files
    pdbqt_files = sorted(glob.glob(str(in_dir / "*.pdbqt")))
    dlg_files   = sorted(glob.glob(str(in_dir / "*.dlg")))
    files = [(f, False) for f in pdbqt_files] + [(f, True) for f in dlg_files]
    if not files:
        sys.exit(f"No .pdbqt or .dlg files found in {in_dir}")

    # Optional SMILES registry
    smiles_map: Dict[str, str] = {}
    if args.smiles_csv:
        df_sm = pd.read_csv(args.smiles_csv)
        cols = {c.lower(): c for c in df_sm.columns}
        if "cid" not in cols or "smiles" not in cols:
            raise ValueError("smiles_csv must contain columns: CID,smiles")
        df_sm = df_sm.rename(columns={cols["cid"]: "CID", cols["smiles"]: "smiles"})
        for _, r in df_sm.iterrows():
            smiles_map[str(r["CID"])] = str(r["smiles"])

    # Writers
    hits_writer     = Chem.SDWriter(str(out_dir / "hits.sdf"))
    top_hits_writer = Chem.SDWriter(str(out_dir / "top_pose_hits.sdf"))

    hits_rows = []
    top_rows  = []

    for fpath, is_dlg in files:
        # CID from filename stem (trim "_dock")
        stem = Path(fpath).stem
        CID = stem_to_cid(stem)

        try:
            mol, scores, pose_idxs = read_pdbqt_with_meeko(fpath, max_poses=args.max_poses, is_dlg=is_dlg)
        except Exception as e:
            print(f"[ERROR] {CID}: {e}", file=sys.stderr)
            continue

        nconf = mol.GetNumConformers()
        if nconf == 0:
            print(f"[WARN] No conformers in {CID}, skipping.")
            continue
        # Compute RMSD matrix
        rms = rmsd_matrix(mol)

        # Write optional RMSD matrix
        if args.write_rmsd_matrices:
            pd.DataFrame(rms, index=pose_idxs, columns=pose_idxs).to_csv(out_dir / f"{CID}_pairwise_rmsd.csv", index=True)

        # Complete-linkage converged cluster (< threshold)
        cluster_idx = greedy_complete_linkage_cluster(rms, scores, args.rmsd_thresh)

        # Hit criteria (≥ hit_min_poses in the converged cluster)
        is_hit = len(cluster_idx) >= args.hit_min_poses
        if not is_hit:
            continue  # Only hits are written

        # Choose top pose within the cluster: lowest docking score
        cluster_sorted = sorted([(i, scores[i]) for i in cluster_idx],
                                key=lambda t: (t[1] if not math.isnan(t[1]) else float("inf")))
        best_idx = cluster_sorted[0][0] if cluster_sorted else None

        # Base properties (once per CID)
        props_base = rdkit_physchem_props(mol)
        smiles = smiles_map.get(CID, mol_to_smiles(mol))

        # ---- HITS: write ONLY viable poses (cluster members) ----
        for i in cluster_idx:
            pose_no = pose_idxs[i]
            sc = scores[i]
            rmsd_to_best = (rms[best_idx, i] if best_idx is not None else float("nan"))
            row = {
                "CID": CID,
                "pose_index": pose_no,
                "docking_score": sc,
                "rmsd_to_best_selected": rmsd_to_best,
                "MWT": props_base["MWT"],
                "TPSA": props_base["TPSA"],
                "logP": props_base["logP"],
                "HBD": props_base["HBD"],
                "HBA": props_base["HBA"],
                "SMILES": smiles,
            }
            hits_rows.append(row)

            m1 = clone_with_conformer(mol, i)
            m1.SetProp("_Name", f"{CID}_pose{pose_no}")
            for k, v in row.items():
                m1.SetProp(k, "" if v is None else str(v))
            hits_writer.write(m1)

        # ---- TOP POSE SUMMARY (one per CID) ----
        stats = cluster_rmsd_stats(rms, cluster_idx, best_idx)
        top_row = {
            "CID": CID,
            "pose_index": pose_idxs[best_idx],
            "docking_score": scores[best_idx],
            "num_viable_poses": len(cluster_idx),
            "rmsd_pairwise_std": stats["rmsd_pairwise_std"],
            "rmsd_to_best_std": stats["rmsd_to_best_std"],
            "MWT": props_base["MWT"],
            "TPSA": props_base["TPSA"],
            "logP": props_base["logP"],
            "HBD": props_base["HBD"],
            "HBA": props_base["HBA"],
            "SMILES": smiles,
        }
        top_rows.append(top_row)

        top_mol = clone_with_conformer(mol, best_idx)
        top_mol.SetProp("_Name", f"{CID}_pose{pose_idxs[best_idx]}")
        for k, v in top_row.items():
            top_mol.SetProp(k, "" if v is None else str(v))
        top_hits_writer.write(top_mol)

    # Close writers
    hits_writer.close()
    top_hits_writer.close()

    # CSVs
    df_hits = pd.DataFrame(hits_rows)
    df_top  = pd.DataFrame(top_rows)

    if not df_hits.empty:
        df_hits = df_hits.sort_values(by=["docking_score", "CID", "pose_index"], ascending=[True, True, True])
    if not df_top.empty:
        df_top = df_top.sort_values(by=["docking_score", "CID"], ascending=[True, True])

    df_hits.to_csv(out_dir / "hits.csv", index=False)
    df_top.to_csv(out_dir / "top_pose_hits.csv", index=False)

    print(
        f"\nDone.\nOutputs in: {out_dir.resolve()}\n"
        f"- hits.sdf / hits.csv (only hits; only viable poses; includes SMILES)\n"
        f"- top_pose_hits.sdf / top_pose_hits.csv (one top pose per hit; includes SMILES; RMSD std metrics)\n"
    )


if __name__ == "__main__":
    main()
