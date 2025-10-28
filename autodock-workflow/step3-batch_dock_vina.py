#!/usr/bin/env python3
"""
 Step 3 :Batch docking with AutoDock Vina Python API.        
 step3-batch_dock_vina.py

 Summary:
 Automates large-scale ligand docking using the AutoDock Vina Python API.
 The script distributes docking jobs across multiple CPU processes for speed,
 reusing precomputed receptor maps within each worker to minimize overhead.
 It reads a defined search box (center/size) from a text file, performs
 ligand minimization and docking, and writes output poses and summary CSVs
 containing best docking energies and job status.

 Key Features:
 - Parallel docking using ProcessPoolExecutor
 - Caching of receptor maps per worker for performance
 - Optional minimization of ligands before docking
 - Automatic resume (skips completed ligands unless rescoring is requested)
 - Live progress tracking with tqdm

 Usage:
   python step3-batch_dock_vina.py \
       --receptor receptor.pdbqt \
       --box docking_box.txt \
       --ligands "ligands/*.pdbqt" \
       --outdir docked_out_api \
       --exhaustiveness 8 --n-poses 5 --processes 8


 Usage help :
      
      python step3-batch_dock_vina.py --help

 Output:
   - Docked ligand files: <outdir>/<ligand>_dock.pdbqt
   - CSV summaries: _manifest_api.csv (all jobs), _scores_api.csv (sorted energies)

"""


import os
import glob
import time
import math
import argparse
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from vina import Vina

# progress bar (optional)
try:
    from tqdm.auto import tqdm
    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False

# -------------------- utilities --------------------

def read_box_txt(path: str):
    """Parse center/size from a simple .txt box file."""
    cx=cy=cz=sx=sy=sz=None
    with open(path) as fh:
        for line in fh:
            t = line.strip().replace(",", " ").replace("=", " ")
            if not t:
                continue
            p  = t.split()
            lo = [w.lower() for w in p]
            if "center_x" in lo: cx = float(p[-1])
            elif "center_y" in lo: cy = float(p[-1])
            elif "center_z" in lo: cz = float(p[-1])
            elif "size_x"   in lo: sx = float(p[-1])
            elif "size_y"   in lo: sy = float(p[-1])
            elif "size_z"   in lo: sz = float(p[-1])
            elif lo[0] == "center" and len(p) >= 4: cx, cy, cz = map(float, p[-3:])
            elif lo[0] == "size"   and len(p) >= 4: sx, sy, sz = map(float, p[-3:])
    assert None not in (cx,cy,cz,sx,sy,sz), f"Bad box file: {path}"
    return (cx,cy,cz), (sx,sy,sz)

def get_best_energy(v: Vina, n_poses: int):
    """
    Version-proof extraction of top docking score.
    Vina.energies() may return (energies, rmsd_lb, rmsd_ub) or (energies, rmsd).
    """
    res = v.energies(n_poses=n_poses)
    energies = res[0] if isinstance(res, (list, tuple)) and len(res) >= 1 else []
    try:
        return float(energies[0][0]) if len(energies) else float("nan")
    except Exception:
        return float("nan")

# -------------------- per-process cache --------------------

_global_vina = None
_global_cfg  = None   # (receptor, center, size, scoring, seed)

def dock_worker(lig_path: str,
                receptor: str,
                center: tuple,
                size: tuple,
                scoring: str,
                seed: int,
                exhaustiveness: int,
                n_poses: int,
                outdir: str,
                rescore_existing: bool = False,
                mute_output: bool = True):
    """
    Worker function that lives in a separate process. Reuses a cached Vina object
    (with precomputed maps) inside the process for speed.
    """
    # Avoid oversubscription (important when running many processes)
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    global _global_vina, _global_cfg
    name = os.path.splitext(os.path.basename(lig_path))[0]
    out_pdbqt = os.path.join(outdir, f"{name}_dock.pdbqt")

    try:
        # Prepare/reuse Vina maps per worker
        cfg = (receptor, center, size, scoring, seed)
        if _global_vina is None or _global_cfg != cfg:
            v = Vina(sf_name=scoring, seed=seed)
            v.set_receptor(receptor)                  # positional arg = version-safe
            v.compute_vina_maps(center=center, box_size=size)
            _global_vina = v
            _global_cfg  = cfg
        v = _global_vina


        # Fresh docking
        v.set_ligand_from_file(lig_path)
        # --- Optional pre-docking local minimization ---
        try:
            energy_minimized = v.optimize()
            min_path = out_pdbqt.replace("_dock.pdbqt", "_minimized.pdbqt")
            v.write_pose(min_path, overwrite=True)
            print(f"[min] {name}: minimized energy {energy_minimized[0]:.2f} kcal/mol")
        except Exception as e:
            print(f"[warn] Minimization failed for {name}: {e}")

        # --- Docking step ---

        if mute_output:
            with open(os.devnull, "w") as _null, contextlib.redirect_stdout(_null), contextlib.redirect_stderr(_null):
                v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        else:
            v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)
        best = get_best_energy(v, n_poses=n_poses)

        return {"CID": name, "status": "ok", "out_pdbqt": out_pdbqt,
                "best_kcal_mol": best, "error": ""}

    except Exception as e:
        return {"CID": name, "status": "error", "out_pdbqt": "",
                "best_kcal_mol": float("nan"), "error": repr(e)[:400]}

# -------------------- main entry --------------------

def main():
    ap = argparse.ArgumentParser(description="Batch docking with AutoDock Vina Python API.")
    ap.add_argument("--receptor", required=True, help="Receptor PDBQT path")
    ap.add_argument("--box",      required=True, help="Box text file (center/size)")
    ap.add_argument("--ligands",  required=True, help="Glob for ligand PDBQT files, e.g. 'pdbqt_ligands_filtered/*.pdbqt'")
    ap.add_argument("--outdir",   default="docked_out_api", help="Output directory")
    ap.add_argument("--scoring",  default="vina", choices=["vina","ad4"], help="Scoring function")
    ap.add_argument("--exhaustiveness", type=int, default=8, help="Vina exhaustiveness")
    ap.add_argument("--n-poses",  type=int, default=5, help="Number of poses to save")
    ap.add_argument("--seed",     type=int, default=42, help="Random seed")
    ap.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 4) - 1),
                    help="Number of worker processes")
    ap.add_argument("--limit",    type=int, default=None, help="Limit to first N ligands for dry run")
    ap.add_argument("--rescore-existing", action="store_true", help="Rescore existing output files instead of skipping")
    ap.add_argument("--unmute", action="store_true", help="Show Vina progress bars from workers")
    ap.add_argument("--run-tag", default=None, help="Optional tag appended to output directory (e.g., run1).")
    ap.add_argument("--clean-run", action="store_true", help="If set, create a new timestamped output directory per run.")
    ap.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar (fallback to prints).")
    args = ap.parse_args()

    # Setup
    os.makedirs(args.outdir, exist_ok=True)
    center, size = read_box_txt(args.box)
    ligs_all = sorted(glob.glob(args.ligands))
    if not ligs_all:
        raise SystemExit(f"No ligands matched: {args.ligands}")
    ligs = ligs_all[:args.limit] if args.limit else ligs_all

    # Cap process count sensibly
    processes = max(1, min(args.processes, (os.cpu_count() or 4)))
    # Throughput generally peaks around physical core count with OMP_NUM_THREADS=1
    print(f"Docking {len(ligs)} ligands (of {len(ligs_all)} total). Using {processes} processes.")

    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=processes) as ex:
        futs = [ex.submit(
            dock_worker,
            p, args.receptor, center, size, args.scoring, args.seed,
            args.exhaustiveness, args.n_poses, args.outdir,
            args.rescore_existing, not args.unmute
        ) for p in ligs]

        use_bar = _HAS_TQDM and (not args.no_progress)
        pbar = None
        if use_bar:
            pbar = tqdm(total=len(ligs), unit="lig", dynamic_ncols=True, leave=True)
            pbar.set_description("Docking")

        # live progress
        for i, f in enumerate(as_completed(futs), 1):
            res = f.result()
            rows.append(res)

            # counts
            oks   = sum(r["status"] == "ok"      for r in rows)
            skip  = sum(r["status"] == "skipped" for r in rows)
            resc  = sum(r["status"] == "rescore" for r in rows)
            errs  = sum(r["status"] == "error"   for r in rows)

            if pbar:
                pbar.update(1)
                # avoid super noisy postfix updates
                if (i % max(5, math.ceil(len(ligs)/50))) == 0 or i == len(ligs):
                    pbar.set_postfix(OK=oks, Resc=resc, Skip=skip, Err=errs)
            else:
                if (i % max(5, math.ceil(len(ligs)/20))) == 0 or i == len(ligs):
                    print(f"[pool] {i}/{len(ligs)}  OK {oks}  Resc {resc}  Skip {skip}  Err {errs}", flush=True)

        if pbar:
            pbar.close()

    dt = time.time() - t0
    print(f"Done {len(rows)} ligands in {dt:.1f}s  (~{len(rows)/max(dt,1):.2f} lig/s)")

    # Manifests
    manifest = pd.DataFrame(rows, columns=["CID","status","out_pdbqt","best_kcal_mol","error"])
    manifest.to_csv(os.path.join(args.outdir, "_manifest_api.csv"), index=False)

    scores = (manifest
              .dropna(subset=["best_kcal_mol"])
              .sort_values("best_kcal_mol"))
    scores.to_csv(os.path.join(args.outdir, "_scores_api.csv"), index=False)

    # Friendly summary
    print("Top 10 by best_kcal_mol:")
    # print(scores.head(10).to_string(index=False))

if __name__ == "__main__":
    main()

