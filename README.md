# CDK2–Staurosporine: Docking, MD and MM-GBSA Workflow

End-to-end structure-based workflow on a well-studied kinase target: Vina docking, GROMACS all-atom MD, and MM-GBSA binding energy with a two-window consistency check.

**Status:** pipeline validated on the native ligand (10 ns test). Other ligands, a longer run and a ligand ranking vs. experiment are **not done yet** (see Limitations).

## System

| Item | Value |
|---|---|
| Target | CDK2, PDB 1AQ1, chain A |
| Ligand | Staurosporine (PDB code STU; residue name `UNL` in the simulation) |
| Simulated system | 50,410 atoms; 46 Na⁺ and 50 Cl⁻ (neutralized, ~0.15 M salt) |

## Methods

> Check every line in this section against the scripts you actually ran, and fix any differences before publishing.

1. **Preparation:** protein chain A only, waters and alternate locations removed; native ligand extracted for box placement and redocking.
2. **Docking:** AutoDock Vina, box centred on the crystal ligand. Redocking RMSD: **[fill in from `dock_results.csv`]** Å.
3. **Ligand parameters:** GAFF2 with AM1-BCC charges (ACPYPE/AmberTools), net charge 0 (set from the input SMILES).
4. **Protein force field:** AMBER99SB-ILDN, TIP3P water, dodecahedral box (1.0 nm margin).
5. **Equilibration:** steepest-descent minimization (converged in 1062 steps), NVT 100 ps, NPT 200 ps with position restraints, 300 K.
6. **Production:** 10 ns, 2 fs timestep, one run, 1001 frames.
7. **Analysis:** GROMACS (`rms`, `rmsf`, `gyrate`) plus MDAnalysis for contacts and H-bonds (`analyze.py`).
8. **MM-GBSA:** gmx_MMPBSA, GB model igb=5, 0.15 M salt, no entropy term. Two non-overlapping windows (2–6 ns and 6–10 ns), 81 frames each.

**Software:** GROMACS 2026.3 (conda-forge), gmx_MMPBSA 1.6.5, MDAnalysis 2.9, Python 3.10. See `environment.yml`.

**Notes for reproducibility:** `-maxwarn 1` was used at the ion-placement `grompp` step (net-charge/excluded-atom warning before minimization). gmx_MMPBSA reported 37 invalid dihedral periodicities = 0 in the converted topology and set them to 1.

## Results (native ligand, 10 ns)

**Equilibration**

| Check | Value |
|---|---|
| NPT temperature | 299.9 K |
| NPT density | 1019.5 kg/m³ |
| Production potential energy | −677,470 kJ/mol (mean) |

**Stability** (`results/analysis.png`)

| Metric | Value |
|---|---|
| Protein backbone RMSD | 2.02 ± 0.26 Å (second-half drift +0.002 Å/ns) |
| Ligand RMSD | 1.33 ± 0.43 Å (second-half drift +0.09 Å/ns) |
| Radius of gyration | 20.41 ± 0.12 Å |
| Min ligand–protein distance | mean 2.81 Å, max 3.17 Å |
| Residues within 4 Å of ligand | 15.9 per frame (mean) |

The ligand stayed in the pocket throughout. It had one transient excursion near 8 ns (ligand RMSD up to ~3.3 Å) that returned to baseline.

**Contacts and hydrogen bonds** (3.0 Å, 150°)

| Interaction | Occupancy |
|---|---|
| Ligand N1 → Glu81 O | 76.0% |
| Leu83 N → ligand O1 | 54.1% |
| Ligand N3 → Gln131 O | 9.7% |

Residues within 4 Å for ≥90% of frames: Glu81, Leu83, Gln131, Asp145, Phe82, Ile10, Ala31, Leu134, Lys33, Phe80.

**MM-GBSA** (kcal/mol, mean; SD across frames in brackets)

| Window | ΔTOTAL | ΔVDWAALS | ΔEEL | ΔEGB | ΔESURF |
|---|---|---|---|---|---|
| 2–6 ns | −48.76 (3.02) | −55.35 | −29.93 | +42.92 | −6.40 |
| 6–10 ns | −47.38 (3.22) | [fill in] | [fill in] | [fill in] | [fill in] |

Both windows agree to within about 1.4 kcal/mol. In window 1 the van der Waals term is the largest favourable contribution; the favourable gas-phase electrostatics (−29.9) are largely cancelled by polar desolvation (+42.9).

## Limitations

- Single 10 ns trajectory, one ligand: no replicates and no ranking against experimental IC50 values.
- MM-GBSA here has no entropy term, so ΔTOTAL is **not** an absolute binding free energy and should not be compared directly to experiment.
- SD values are across correlated frames, so they understate the true uncertainty.
- Ligand protonation was fixed from the input SMILES.

## Planned work

- Dock and simulate additional CDK2 inhibitors with literature IC50 values and test whether MM-GBSA ranking matches experiment.
- Longer native run and replicate simulations.

## Repository layout

```
scripts/    dock.py, make_lig.py, run_md.sh, analyze.py, collect.py
mdp/        GROMACS parameter files
results/    analysis.png, summary.txt, results_win1.dat, results_win2.dat, dock_results.csv
environment.yml
```

Large files (trajectories, checkpoints, `.tpr`) are not tracked; see `.gitignore`.
