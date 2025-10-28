#!/usr/bin/env python3

"""
 Author : Payal Chatterjee

 Script Name: step2-prep-ligs.py

 Summary:
 Prepares ligand libraries for virtual screening by fetching, cleaning, and filtering compounds 
 similar to a co-crystal ligand from RCSB. 
 Integrates PubChem and Enamine sources, cleans molecules using molscrub at pH 7.4, 
 applies physicochemical filters (MW, TPSA, logP, rotatable bonds, halogens) 
 and structural filters (BRENK, NIH, and custom alert patterns), computes ECFP4 similarities to the co-ligand, 
 removes near-duplicates, and exports the top curated set as <protein>_ligands_<n>.sdf and Meeko-ready PDBQT files for downstream docking.
 
 Workflow / Key Features:
 - Fetch co-crystal ligand and canonical SMILES from RCSB
 - Identify structurally similar analogs from PubChem (Tanimoto ≥ threshold)
 - Clean all molecules using molscrub at pH 7.4
 - Apply physicochemical and structural filters (BRENK, NIH, and custom alerts)
 - Process Enamine library in parallel batches with molscrub filtering
    - could alternatively be done using meeko + openbabel(for protonation states) + RDKit : would be faster but more no. of steps/sanitation bugs
 - Compute ECFP4 similarities to the co-ligand and remove near-duplicates
 - Rank and retain top N diverse ligands for screening
 - Export cleaned library as <protein>_ligands_<n>.sdf
 - Generate Meeko-ready PDBQT files with manifest for AutoDock/DiffDock workflows

 Usage:

python step2-prep-ligs.py \
       --pdb-id 1H1Q --lig-id 2A6 --chain-id A \
       --pubchem-threshold 85 --pubchem-max 50 \
       --enamine-sdf Enamine_Hinge_Binders.sdf \
       --outdir results --n-jobs 8 --drop-alerts

 Output:
   - Reference ligand files: <lig-id>_ref_as_modeled.mol / _ref_protonated.mol / .pdb
   - Filtered ligand library: <pdb-id>_ligands_<n>.sdf
   - Prepared docking inputs: pdbqt_scrubbed/ (PDBQT + _manifest.csv)
   - Optional figures: figs/ (similarity histograms and ECDFs)

 Notes:
   The final SDF naming convention follows <protein>_ligands_<count>.sdf,
   based on user-provided --pdb-id and number of ligands retained.

"""

from __future__ import annotations

import os
import sys
import io
import gzip
import contextlib
import argparse
import warnings

# Use non-interactive Matplotlib backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# External deps
import requests
import pubchempy as pcp
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# RDKit
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors, Crippen, FilterCatalog, PandasTools, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: F401 (import preserved)
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.Chem.SaltRemover import SaltRemover         # noqa: F401 (import preserved)

# Silence RDKit logs
RDLogger.DisableLog("rdApp.*")

# Warnings: prody deprecation filter (preserved behavior)
os.environ["PYTHONWARNINGS"] = "ignore:pkg_resources is deprecated as an API.*:UserWarning:prody.utilities.misctools"
warnings.filterwarnings(
    "ignore",
    message=r"pkg_resources is deprecated as an API.*",
    category=UserWarning,
    module=r"prody\.utilities\.misctools"
)

# May be optional in some environments; keep import as in original flow
try:
    import prody  # noqa: F401
except Exception:
    pass

# --- Globals from the notebook (preserved) ---

def fetch_ligand_instance_sdf(pdb_id: str, lig_id: str, chain_id: str | None = None):
    """
    Download ligand instance(s) from RCSB ModelServer as SDF with correct bonds/charges.
    Returns a list of RDKit Mol objects (each with a 3D conformer from the crystal).
    """
    base = f"https://models.rcsb.org/v1/{pdb_id}/ligand?label_comp_id={lig_id}&encoding=sdf"
    url = base + (f"&auth_asym_id={chain_id}" if chain_id else "")
    r = requests.get(url)
    r.raise_for_status()

    # RDKit SDF supplier on bytes
    supp = Chem.ForwardSDMolSupplier(io.BytesIO(r.content), removeHs=False, sanitize=True)
    mols = [m for m in supp if m is not None]
    if not mols:
        raise ValueError(f"No ligand molecules returned for {pdb_id}:{lig_id} (chain={chain_id})")

    # Optional: report instance identifiers if present
    for i, m in enumerate(mols):
        name = m.GetProp("_Name") if m.HasProp("_Name") else f"{lig_id}"
        print(f"[INFO] Instance {i}: name={name}, atoms={m.GetNumAtoms()}, confs={m.GetNumConformers()}")
    return mols


def add_hydrogens_preserve_coords(m: Chem.Mol) -> Chem.Mol:
    """
    Adds explicit hydrogens using existing 3D coordinates; heavy atoms are unchanged.
    """
    if m.GetNumConformers() == 0:
        raise ValueError("Ligand has no 3D conformer; cannot place hydrogens with coordinates.")
    mH = Chem.AddHs(m, addCoords=True)
    return mH


def get_smiles_from_rcsb_graphql(ligand_ids):
    """Return a {ligand_id: canonical_SMILES} dict from RCSB GraphQL API."""
    base_url = "https://data.rcsb.org/graphql"
    query_fmt = "[" + ", ".join(['"' + i + '"' for i in ligand_ids]) + "]"
    query = f"""
    {{
      chem_comps(comp_ids:{query_fmt}) {{
        chem_comp {{ id }}
        rcsb_chem_comp_descriptor {{ SMILES_stereo }}
      }}
    }}
    """
    r = requests.post(base_url, json={"query": query})
    r.raise_for_status()
    data = r.json()["data"]["chem_comps"]
    out = {}
    for d in data:
        cid = d["chem_comp"]["id"]
        desc = d.get("rcsb_chem_comp_descriptor")
        if desc and desc.get("SMILES_stereo"):
            out[cid] = desc["SMILES_stereo"]
    return out


def get_pubchem_similar(smiles, threshold=85, max_results=50):
    """Fetch PubChem compounds similar to input SMILES (Tanimoto >= threshold)."""
    cids = pcp.get_cids(smiles, namespace='smiles', searchtype='similarity',
                        threshold=threshold, listkey_count=max_results)
    comps = pcp.get_compounds(cids[:max_results])
    data = []
    for c in comps:
        data.append({
            "CID": c.cid,
            #"SMILES": c.smiles,
            "SMILES": getattr(c, "canonical_smiles", getattr(c, "isomeric_smiles", None)),
            "IUPACName": c.iupac_name,
            "MolWt": c.molecular_weight,
            "LogP": c.xlogp,
            "RotBonds": c.rotatable_bond_count,
            "TPSA": c.tpsa
        })
    return pd.DataFrame(data)


# -- BRENK + NIH + custom filters (preserved as-is) --
CUSTOM_ALERTS = {
    "Carboxylate": "[CX3](=O)[O- ", "AcidHalide": "C(=O)[Cl,Br,F,I]",
    "AcidAzide": "C(=O)N=[N+]=[N-]",
    "SulfonylHalide": "S(=O)(=O)[Cl,Br,F,I]",
    "SulfonylAzide": "S(=O)(=O)N=[N+]=[N-]",
    "Catechol": "c1ccc(O)c(O)c1",
    "Isocyanate": "N=C=O",
}
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
FILTER = FilterCatalog.FilterCatalog(params)
CUSTOM_PATTS = {n: Chem.MolFromSmarts(s) for n, s in CUSTOM_ALERTS.items()}


def alert_hits(m):
    hits = [e.GetDescription() for e in FILTER.GetMatches(m)]
    for name, patt in CUSTOM_PATTS.items():
        if patt and m.HasSubstructMatch(patt):
            hits.append(name)
    return sorted(set(hits))


def strict_rotb(m): return int(rdMolDescriptors.CalcNumRotatableBonds(m, strict=True))
def halogen_count(m): return sum(a.GetAtomicNum() in (9, 17, 35, 53) for a in m.GetAtoms())


def passes_physchem(m,
                    mw_range=(250, 600),
                    tpsa_range=(40, 140),
                    logp_max=5.0,
                    rotb_max=10,
                    hal_max=2):
    try:
        mw = Descriptors.MolWt(m)
        logp = Descriptors.MolLogP(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        rotb = int(rdMolDescriptors.CalcNumRotatableBonds(m, strict=True))
        halo = sum(a.GetAtomicNum() in (9, 17, 35, 53) for a in m.GetAtoms())
        ok = (
            mw_range[0] <= mw <= mw_range[1]
            and tpsa_range[0] <= tpsa <= tpsa_range[1]
            and logp < logp_max
            and rotb <= rotb_max
            and halo <= hal_max
        )
        return ok, (mw, logp, tpsa, rotb, halo)
    except Exception:
        return False, (None, None, None, None, None)


@contextlib.contextmanager
def hush_all():
    from rdkit import RDLogger, rdBase
    RDLogger.DisableLog('rdApp.*')
    with rdBase.BlockLogs():
        with open(os.devnull, "w") as devnull, \
             contextlib.redirect_stdout(devnull), \
             contextlib.redirect_stderr(devnull):
            yield


# --- molscrub setup (preserved usage) ---
try:
    from molscrub import Scrub
    scrubber = Scrub(ph_low=7.4, ph_high=7.4)
except Exception as e:
    print("[ERROR] molscrub is required for this workflow. Install with `pip install molscrub`.")
    raise
def _clean_one_pubchem_molscrub(row, mw_range, tpsa_range, logp_max, rotb_max, hal_max, pH, drop_alerts):
    cid = row["CID"]; smi = row["SMILES"]
    if not isinstance(smi, str) or not smi:
        return None

    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None

        scrubbed_states = list(scrubber(m))
        if not scrubbed_states:
            return None
        m = scrubbed_states[0]
    except Exception:
        return None

    hits = alert_hits(m)
    if drop_alerts and hits:
        return None

    ok, (mw, logp, tpsa, rotb, halo) = passes_physchem(
        m, mw_range, tpsa_range, logp_max, rotb_max, hal_max
    )
    if not ok:
        return None

    smi_out = Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
    return {
        "CID": cid,
        "SMILES": smi_out,
        "MolWt": mw,
        "logP": logp,
        "TPSA": tpsa,
        "RotBonds": rotb,
        "HalogenCount": halo,
        "AlertHits": ";".join(hits),
        "AlertAny": bool(hits)
    }


def _process_pubchem_chunk_molscrub(rows, mw_range, tpsa_range, logp_max, rotb_max, hal_max, pH, drop_alerts):
    return [
        _clean_one_pubchem_molscrub(r, mw_range, tpsa_range, logp_max, rotb_max, hal_max, pH, drop_alerts)
        for r in rows
    ]


def clean_pubchem_smiles_parallel_molscrub(analogs_df: pd.DataFrame,
                                           mw_range=(250, 600), tpsa_range=(40, 140),
                                           logp_max=5.0, rotb_max=10, hal_max=2,
                                           pH=7.4, n_jobs=None, chunksize=256, drop_alerts=False):
    """Parallel cleaning of PubChem molecules using Molscrub Python API."""
    n_jobs = n_jobs or os.cpu_count() or 4
    rows = analogs_df[["CID", "SMILES"]].dropna().to_dict("records")
    out = []
    with hush_all():
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futs = []
            for i in range(0, len(rows), chunksize):
                futs.append(ex.submit(_process_pubchem_chunk_molscrub,
                                      rows[i:i+chunksize],
                                      mw_range, tpsa_range, logp_max, rotb_max, hal_max, pH, drop_alerts))
            for f in tqdm(as_completed(futs), total=len(futs), desc="PubChem", unit="chunk"):
                out.extend([r for r in f.result() if r is not None])
    df = pd.DataFrame(out, columns=[
        "CID", "SMILES", "MolWt", "logP", "TPSA", "RotBonds", "HalogenCount", "AlertHits", "AlertAny"
    ])
    if not df.empty:
        df["Source"] = "PubChem"
    print(f"PubChem kept: {len(df)} / {len(analogs_df)}")
    return df


def _yield_sdf_batches(path, batch_size=500, sanitize=True, removeHs=False):
    """
    Yield batches of valid molecule MolBlocks and CIDs from an SDF file.
    MolBlocks (strings) are lightweight and pickle-safe for multiprocessing.
    """
    suppl = Chem.ForwardSDMolSupplier(path, sanitize=False, removeHs=removeHs)
    molblocks, cids = [], []

    for m in suppl:
        if m is None:
            continue

        if sanitize:
            try:
                Chem.SanitizeMol(m, catchErrors=True)
            except Exception:
                continue

        cid = None
        for key in ["CID", "Catalog ID", "_Name", "ID", "Compound_ID"]:
            if m.HasProp(key):
                cid = m.GetProp(key).strip()
                break
        if cid is None:
            cid = f"UNK_{len(molblocks)+1}"

        molblocks.append(Chem.MolToMolBlock(m))
        cids.append(cid)

        if len(molblocks) >= batch_size:
            yield molblocks, cids
            molblocks, cids = [], []

    if molblocks:
        yield molblocks, cids


def _enamine_worker_batch_molscrub(molblocks, cids,
                                   mw_range, tpsa_range,
                                   logp_max, rotb_max, hal_max,
                                   pH, drop_alerts):
    out = []
    bad = propdrop = phfixed = 0

    for mb, cid in zip(molblocks, cids):
        try:
            m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
            if m is None:
                bad += 1
                continue

            Chem.SanitizeMol(m, catchErrors=True)

            scrubbed_states = list(scrubber(m))
            if not scrubbed_states:
                bad += 1
                continue
            m = scrubbed_states[0]
            phfixed += 1
        except Exception:
            bad += 1
            continue

        hits = alert_hits(m)
        if drop_alerts and hits:
            propdrop += 1
            continue

        ok, (mw, logp, tpsa, rotb, halo) = passes_physchem(
            m, mw_range, tpsa_range, logp_max, rotb_max, hal_max
        )
        if not ok:
            propdrop += 1
            continue

        smi = Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
        out.append((cid, smi, mw, logp, tpsa, rotb, halo, ";".join(hits), bool(hits)))

    return out, bad, propdrop, phfixed


def clean_enamine_parallel_molscrub(path,
                                    mw_range=(250, 600), tpsa_range=(40, 140),
                                    logp_max=5.0, rotb_max=10, hal_max=2,
                                    pH=7.4, batch_size=500, n_jobs=None, drop_alerts=False):
    """
    Fully parallel Enamine cleanup using molscrub.
    Reads SDF in small batches and distributes to multiple processes.
    """
    n_jobs = n_jobs or os.cpu_count() or 4
    results = []
    total_bad = total_propdrop = total_phfixed = 0

    with hush_all():
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futs = []
            for molblocks, cids in _yield_sdf_batches(path, batch_size=batch_size):
                futs.append(ex.submit(_enamine_worker_batch_molscrub,
                                      molblocks, cids,
                                      mw_range, tpsa_range,
                                      logp_max, rotb_max, hal_max,
                                      pH, drop_alerts))

            for fut in tqdm(as_completed(futs), total=len(futs), desc="Enamine", unit="batch"):
                out, bad, propdrop, phfixed = fut.result()
                total_bad += bad
                total_propdrop += propdrop
                total_phfixed += phfixed
                results.extend(out)

    df = pd.DataFrame(results, columns=[
        "CID", "SMILES", "MolWt", "logP", "TPSA", "RotBonds", "HalogenCount", "AlertHits", "AlertAny"
    ])
    if not df.empty and df["CID"].isna().any():
        miss = df["CID"].isna()
        df.loc[miss, "CID"] = [f"ENM_{i}" for i in range(1, miss.sum()+1)]
    if not df.empty:
        df["Source"] = "Enamine"

    print(f"Enamine kept: {len(df)} | invalid: {total_bad} | drops: {total_propdrop} | pH-fixed: {total_phfixed}")
    return df


def canonicalize_smiles(smi: str):
    """Return canonical, isomeric SMILES from input string."""
    try:
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(m, isomericSmiles=True, canonical=True) if m else None
    except Exception:
        return None


def prep_with_scrub(smi, lig_id="Lig"):
    """Return cleaned, protonated RDKit Mol (3D) from SMILES via Molscrub."""
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return None
    try:
        states = list(scrubber(mol))
        if not states:
            return None
        m = states[0]
        Chem.SanitizeMol(m, catchErrors=True)
        m.SetProp("_Name", lig_id)
        return m
    except Exception as e:
        print(f"[WARN] Molscrub failed for {lig_id}: {e}")
        return None


def write_pdbqt_from_df(df, mol_col="Mols", id_col="CID",
                        outdir="pdbqt_ready", charge_model="gasteiger"):
    """
    Convert protonated, 3D RDKit Mol objects to PDBQT via Meeko.
    Handles both legacy and new Meeko return types.
    """
    os.makedirs(outdir, exist_ok=True)

    # Meeko imports preserved
    try:
        from meeko import PDBQTWriterLegacy as PDBQTWriter
    except ImportError:
        from meeko import PDBQTWriter

    from meeko import MoleculePreparation

    prep = MoleculePreparation(charge_model=charge_model)
    writer = PDBQTWriter()

    n_ok = n_fail = 0
    rows = []

    for _, row in df.iterrows():
        lig_id = str(row.get(id_col, "Lig"))
        m = row.get(mol_col, None)
        if not isinstance(m, Chem.Mol):
            rows.append({"name": lig_id, "status": "fail", "pdbqt": "", "error": "No Mol object"})
            n_fail += 1
            continue

        try:
            # Ensure explicit Hs
            if not any(a.GetAtomicNum() == 1 for a in m.GetAtoms()):
                m = Chem.AddHs(m, addCoords=True)

            # Ensure conformer
            if m.GetNumConformers() == 0:
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                cid = AllChem.EmbedMolecule(m, params)
                if cid != 0:
                    raise RuntimeError("3D embedding failed")
                AllChem.MMFFOptimizeMolecule(m, maxIters=200)

            # Meeko preparation (returns MoleculeSetup or list thereof)
            setup = prep.prepare(m)
            setup_list = setup if isinstance(setup, (list, tuple)) else [setup]

            # Write one or more PDBQTs
            pdbqt_file = ""
            for i, ms in enumerate(setup_list):
                pdbqt_file = os.path.join(outdir, f"{lig_id}.pdbqt" if i == 0 else f"{lig_id}__alt{i+1}.pdbqt")
                pdbqt_str = writer.write_string(ms)
                if isinstance(pdbqt_str, tuple):
                    pdbqt_str = pdbqt_str[0]
                with open(pdbqt_file, "w") as fh:
                    fh.write(pdbqt_str)
            n_ok += 1
            rows.append({"name": lig_id, "status": "ok", "pdbqt": pdbqt_file, "error": ""})

        except Exception as e:
            n_fail += 1
            rows.append({"name": lig_id, "status": "fail", "pdbqt": "", "error": str(e)})
            print(f"[fail] {lig_id}: {e}")

    manifest = pd.DataFrame(rows)
    manifest.to_csv(os.path.join(outdir, "_manifest.csv"), index=False)
    print(f"[done] Meeko write: OK {n_ok} | Fail {n_fail} → {outdir}")
    return manifest


def parse_args():
    p = argparse.ArgumentParser(
        description="Execute the VS pipeline (notebook -> CLI) without changing computation flow."
    )
    # Step 1 inputs
    p.add_argument("--pdb-id", default="1H1Q", help="PDB ID (default: 1H1Q)")
    p.add_argument("--lig-id", default="2A6", help="Ligand ID (default: 2A6)")
    p.add_argument("--chain-id", default="A", help="Chain ID (default: A, use empty to pull all)")

    # Step 2 PubChem similarity
    p.add_argument("--pubchem-threshold", type=int, default=85, help="PubChem similarity threshold (default: 85)")
    p.add_argument("--pubchem-max", type=int, default=50, help="Max PubChem results (default: 50)")

    # Step 4 filters
    p.add_argument("--mw-min", type=float, default=250)
    p.add_argument("--mw-max", type=float, default=600)
    p.add_argument("--tpsa-min", type=float, default=40)
    p.add_argument("--tpsa-max", type=float, default=140)
    p.add_argument("--logp-max", type=float, default=5.0)
    p.add_argument("--rotb-max", type=int, default=10)
    p.add_argument("--hal-max", type=int, default=2)
    p.add_argument("--drop-alerts", action="store_true", help="Drop molecules with BRENK/NIH/custom alerts")

    # Parallelism / batching
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--pubchem-chunk", type=int, default=256)
    p.add_argument("--enamine-batch-size", type=int, default=1000)

    # Enamine input
    p.add_argument("--enamine-sdf", default="Enamine_Kinase_Hinge_Binders_Library_plated_24000cmpds_20251019.sdf",
                   help="Path to Enamine SDF")

    # Output control
    p.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument("--save-plots", action="store_true", help="Save plots instead of (or in addition to) showing")

    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # --- Step 1: Fetch co-crystal ligand + hydrogens + CCD SMILES (preserved) ---
    pdb_id = args.pdb_id
    lig_id = args.lig_id
    chain_id = args.chain_id if args.chain_id else None

    mols = fetch_ligand_instance_sdf(pdb_id, lig_id, chain_id)
    ref = mols[0]

    Chem.MolToMolFile(ref, os.path.join(args.outdir, f"{lig_id}_ref_as_modeled.mol"))

    refH = add_hydrogens_preserve_coords(ref)
    Chem.MolToMolFile(refH, os.path.join(args.outdir, f"{lig_id}_ref_protonated.mol"))
    Chem.MolToPDBFile(refH, os.path.join(args.outdir, f"{lig_id}_ref_protonated.pdb"))
    print(f"[INFO] Saved {lig_id}_ref_protonated.mol with {refH.GetNumAtoms()} atoms")

    smiles_dict = get_smiles_from_rcsb_graphql([lig_id])

    # display image if possible (preserved intent); also save if requested
    co_ligand_smiles = None
    if lig_id in smiles_dict:
        co_ligand_smiles = smiles_dict[lig_id]
        print(f"RCSB returned canonical SMILES for {lig_id}:")
        print(co_ligand_smiles)

        ref_mol_2d = Chem.MolFromSmiles(co_ligand_smiles)
        if ref_mol_2d:
            img = Draw.MolToImage(ref_mol_2d, size=(250, 250))
            # Try to display if available
            try:
                from IPython.display import display  # noqa
                display(img)
            except Exception:
                pass
            if args.save_plots:
                figs_dir = os.path.join(args.outdir, "figs")
                os.makedirs(figs_dir, exist_ok=True)
                img.save(os.path.join(figs_dir, f"{lig_id}_ccd_2D.png"))
    else:
        print(f"[WARN] No SMILES_stereo descriptor returned for {lig_id}")
        sys.exit(1)

    # --- Step 2: PubChem similar ---
    analogs_df = get_pubchem_similar(co_ligand_smiles,
                                     threshold=args.pubchem_threshold,
                                     max_results=args.pubchem_max)
    print(f"Retrieved {len(analogs_df)} similar compounds from PubChem.")

    ## Visualization grid (preserved generation; save if requested)
    ##mols_grid = [Chem.MolFromSmiles(smi) for smi in analogs_df["SMILES"]]

    #legends = [f"CID {cid}" for cid in analogs_df["CID"]]
    #grid_img = Draw.MolsToGridImage(
    #    mols_grid, molsPerRow=5, subImgSize=(200, 200), legends=legends, useSVG=False
    #)
    #if args.save_plots:
    #    figs_dir = os.path.join(args.outdir, "figs")
    #    os.makedirs(figs_dir, exist_ok=True)
    #    grid_path = os.path.join(figs_dir, "pubchem_similars_grid.png")
    #    grid_img.save(grid_path)

    # --- Step 4: Cleaning (preserved functions and flow) ---
    mw_range = (args.mw_min, args.mw_max)
    tpsa_range = (args.tpsa_min, args.tpsa_max)

    pubchem_clean_df = clean_pubchem_smiles_parallel_molscrub(
        analogs_df,
        mw_range=mw_range, tpsa_range=tpsa_range,
        logp_max=args.logp_max, rotb_max=args.rotb_max, hal_max=args.hal_max,
        n_jobs=args.n_jobs, chunksize=args.pubchem_chunk, drop_alerts=args.drop_alerts
    )

    if not os.path.isfile(args.enamine_sdf):
        print(f"[ERROR] Enamine SDF not found: {args.enamine_sdf}")
        sys.exit(2)

    enamine_clean_df = clean_enamine_parallel_molscrub(
        args.enamine_sdf,
        mw_range=(max(300, args.mw_min), args.mw_max),  # preserve notebook’s stricter MW for Enamine
        tpsa_range=tpsa_range,
        logp_max=args.logp_max, rotb_max=args.rotb_max, hal_max=args.hal_max,
        n_jobs=args.n_jobs, batch_size=args.enamine_batch_size, drop_alerts=args.drop_alerts
    )

    combined = pd.concat([pubchem_clean_df, enamine_clean_df], ignore_index=True)
    combined = combined.dropna(subset=["SMILES"]).drop_duplicates(subset=["SMILES"]).reset_index(drop=True)
    print("Sources\n", combined["Source"].value_counts())
    print("Rows:", len(combined))

    combined_clean = combined[(combined.AlertAny == False)]

    # --- Similarity analysis (preserved constants/flow) ---
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    SMI_COL = "SMILES_CAN" if "SMILES_CAN" in combined_clean.columns else "SMILES"
    ID_COL = "LIG_ID" if "LIG_ID" in combined_clean.columns else "CID"

    if "Source" not in combined_clean.columns:
        combined_clean = combined_clean.copy()
        combined_clean["Source"] = np.where(
            combined_clean[ID_COL].astype(str).str.startswith(("ENM_", "ENAM")),
            "Enamine", "PubChem"
        )

    co_mol = Chem.MolFromSmiles(co_ligand_smiles)
    assert co_mol is not None, "co_ligand_smiles failed to parse."
    co_fp = morgan_gen.GetFingerprint(co_mol)

    def mfp4(smi: str):
        m = Chem.MolFromSmiles(smi)
        return morgan_gen.GetFingerprint(m) if m else None

    def rot_bonds_strict(smi: str):
        m = Chem.MolFromSmiles(smi)
        return int(rdMolDescriptors.CalcNumRotatableBonds(m, strict=True)) if m else np.nan

    df = combined_clean.copy()

    if "FP" not in df.columns:
        df["FP"] = df[SMI_COL].map(mfp4)
    df = df.dropna(subset=["FP"]).reset_index(drop=True)

    if "RotBonds" not in df.columns or df["RotBonds"].isna().any():
        df["RotBonds"] = df[SMI_COL].map(rot_bonds_strict)

    fps = list(df["FP"])
    df["Sim_to_CoLig"] = BulkTanimotoSimilarity(co_fp, fps)

    DROP_NEAR_DUP = True
    SIM_DROP_THR = 0.99
    if DROP_NEAR_DUP:
        before = len(df)
        df = df[df["Sim_to_CoLig"] <= SIM_DROP_THR].reset_index(drop=True)
        print(f"Removed {before - len(df)} ligands with similarity > {SIM_DROP_THR:.2f} to the co-ligand.")

    CLOSE_THR = 0.80
    summary = (
        df.assign(Close=df["Sim_to_CoLig"] >= CLOSE_THR)
          .groupby(["Source", "Close"])[ID_COL]
          .count()
          .unstack(fill_value=0)
          .rename(columns={False: "<= threshold", True: ">= threshold"})
    )
    print(f"\nCounts per source with Sim_to_CoLig >= {CLOSE_THR:.2f}:")
    print(summary)

    # --- Plots (preserved logic; saved if requested) ---
    df_plot = df.copy()
    pc = df_plot.loc[df_plot["Source"] == "PubChem", "Sim_to_CoLig"].dropna()
    en = df_plot.loc[df_plot["Source"] == "Enamine", "Sim_to_CoLig"].dropna()

    bins = np.linspace(0, 1, 41)

    # Similarity histogram
    fig, ax = plt.subplots()
    if len(en): ax.hist(en, bins=bins, histtype="step", density=True, linewidth=2.2, label=f"Enamine (n={len(en)})")
    if len(pc): ax.hist(pc, bins=bins, histtype="step", density=True, linewidth=2.2, label=f"PubChem (n={len(pc)})")
    ax.axvline(CLOSE_THR, ls="--", lw=1.5, label=f"Close thr = {CLOSE_THR:.2f}")
    ax.set_xlabel("Tanimoto (ECFP4) to co-ligand"); ax.set_ylabel("Density")
    ax.legend()
    if args.save_plots:
        figs_dir = os.path.join(args.outdir, "figs")
        os.makedirs(figs_dir, exist_ok=True)
        fig.savefig(os.path.join(figs_dir, "similarity_hist.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ECDF
    def ecdf(x):
        x = np.sort(x)
        y = np.arange(1, len(x)+1) / len(x)
        return x, y

    fig, ax = plt.subplots()
    if len(en):
        xe, ye = ecdf(en.values)
        ax.plot(xe, ye, lw=2, label=f"Enamine (n={len(en)})")
    if len(pc):
        xp, yp = ecdf(pc.values)
        ax.plot(xp, yp, lw=2, label=f"PubChem (n={len(pc)})")
    ax.axvline(CLOSE_THR, ls="--", lw=1.5)
    ax.set_xlabel("Tanimoto (ECFP4) to co-ligand"); ax.set_ylabel("ECDF")
    ax.legend()
    if args.save_plots:
        figs_dir = os.path.join(args.outdir, "figs")
        os.makedirs(figs_dir, exist_ok=True)
        fig.savefig(os.path.join(figs_dir, "similarity_ecdf.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Rotatable-bonds histogram
    pc_rb = df_plot.loc[df_plot["Source"] == "PubChem", "RotBonds"].dropna()
    en_rb = df_plot.loc[df_plot["Source"] == "Enamine", "RotBonds"].dropna()
    rb_bins = np.arange(0, int(df_plot["RotBonds"].max()) + 2)

    fig, ax = plt.subplots()
    if len(en_rb): ax.hist(en_rb, bins=rb_bins, histtype="step", density=True, linewidth=2.2, label="Enamine")
    if len(pc_rb): ax.hist(pc_rb, bins=rb_bins, histtype="step", density=True, linewidth=2.2, label="PubChem")
    ax.axvline(args.rotb_max, ls="--", lw=1.5, label=f"RotBonds cutoff = {args.rotb_max}")
    ax.set_xlabel("# Rotatable Bonds (strict)"); ax.set_ylabel("Density")
    ax.legend()
    if args.save_plots:
        figs_dir = os.path.join(args.outdir, "figs")
        os.makedirs(figs_dir, exist_ok=True)
        fig.savefig(os.path.join(figs_dir, "rotb_hist.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    filtered_df = df.copy()

    # --- Keep notebook’s notes and “top N” selection behavior (preserved) ---
    CLOSE_THR = 0.80
    SIM_DROP_THR = 0.99

    df_before_drop = df.copy()
    n_close_before = (df_before_drop["Sim_to_CoLig"] >= CLOSE_THR).sum()
    n_dropped = (df_before_drop["Sim_to_CoLig"] > SIM_DROP_THR).sum()
    print(f"Close (≥ {CLOSE_THR:.2f}) BEFORE drop: {n_close_before}")
    print(f"Dropped (> {SIM_DROP_THR:.2f}) BEFORE drop: {n_dropped}")

    df_after_drop = df_before_drop[df_before_drop["Sim_to_CoLig"] <= SIM_DROP_THR].reset_index(drop=True)

    TOPN = 200
    df_top = df_before_drop.sort_values("Sim_to_CoLig", ascending=False).head(TOPN)

    # Co-ligand prep
    co_can = canonicalize_smiles(co_ligand_smiles)
    if not co_can:
        raise ValueError("co_ligand_smiles failed to parse.")

    co_mol = prep_with_scrub(co_can, "COLIG_2A6")
    if co_mol is None:
        raise ValueError("Molscrub failed to prepare co-ligand.")

    mw = Descriptors.MolWt(co_mol)
    logp = Crippen.MolLogP(co_mol)
    tpsa = rdMolDescriptors.CalcTPSA(co_mol)
    rotb = int(rdMolDescriptors.CalcNumRotatableBonds(co_mol, strict=True))
    halogens = sum(a.GetAtomicNum() in (9, 17, 35, 53) for a in co_mol.GetAtoms())
    alerts = alert_hits(co_mol)

    # Ensure columns
    df_base = df_top.copy()
    SMI_COL = "SMILES_CAN" if "SMILES_CAN" in df_base.columns else "SMILES"
    ID_COL = "LIG_ID" if "LIG_ID" in df_base.columns else "CID"

    co_row = {
        ID_COL: "COLIG_2A6",
        SMI_COL: co_can,
        "SMILES_CAN": co_can,
        "Source": "CoLigand",
        "IsCoLigand": True,
        "Sim_to_CoLig": 1.0,
        "MolWt": mw,
        "logP": logp,
        "TPSA": tpsa,
        "RotBonds": rotb,
        "HalogenCount": halogens,
        "AlertHits": ";".join(sorted(set(alerts))),
        "AlertAny": bool(alerts),
        "SMILES": co_ligand_smiles,
    }

    if "SMILES_CAN" not in df_base.columns:
        df_base["SMILES_CAN"] = df_base[SMI_COL].map(canonicalize_smiles)

    df_with_co = pd.concat([pd.DataFrame([co_row]), df_base], ignore_index=True)
    df_with_co = df_with_co.drop_duplicates(subset=["SMILES_CAN"], keep="first").reset_index(drop=True)

    print(f"Added co-ligand. Final size: {len(df_with_co)}")

    # Write SDF (preserved behavior)
    df_with_co['Mols'] = [Chem.MolFromSmiles(smic) for smic in df_with_co.SMILES_CAN]
    df_with_co.drop(columns=['IsCoLigand', 'AlertHits', 'AlertAny'], inplace=True)
    #sdf_out = os.path.join(args.outdir, 'cdk2_ligands_top200.sdf')
    sdf_out = os.path.join(args.outdir, f"{pdb_id}_ligands_prepped.sdf")
    PandasTools.WriteSDF(df_with_co, sdf_out, molColName='Mols', idName='CID', properties=df_with_co.columns)
    print(f"[INFO] Wrote SDF: {sdf_out}")

    # PDBQT conversion (preserved)
    pdbqt_dir = os.path.join(args.outdir, "pdbqt_scrubbed")
    manifest = write_pdbqt_from_df(df_with_co, mol_col="Mols", id_col="CID", outdir=pdbqt_dir)

    print("[DONE]")


if __name__ == "__main__":
    main()

