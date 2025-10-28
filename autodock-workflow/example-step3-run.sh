python step3-batch_dock_vina.py --receptor test/1H1Q_receptorH.pdbqt --box test/1H1Q_receptorH.box.txt --ligands test/pdbqt_scrubbed/ --outdir test/outputs/step3-outputs --exhaustiveness 32 --n-poses 5 --processes 8  --limit 20 --clean-run   Docking 1 ligands (of 1 total). Using 8 processes.
Docking:   0%|                                                                                                                                | 0/1 [00:00<?, ?lig/s]Computing Vina grid ... done.
Docking: 100%|███████████████████████████████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  1.41lig/s, Err=1, OK=0, Resc=0, Skip=0]
Done 1 ligands in 0.7s  (~1.00 lig/s)
Top 10 by best_kcal_mol:

