"""
Usage: python analyze.py md/native [LIGAND_RESNAME]
Needs in the folder: md.gro, md_fit.xtc, rmsd_bb.xvg, rmsd_lig.xvg, rmsf.xvg, rg.xvg
"""

import sys, traceback
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import capped_distance

d = sys.argv[1]
xvg = lambda f: np.loadtxt(f"{d}/{f}", comments=("#", "@"))
out = []  # lines saved to summary.txt


def say(s):
    print(s)
    out.append(s)


# ---------- 1. numbers from GROMACS files ----------
rb, rl, rf, rg = xvg("rmsd_bb.xvg"), xvg("rmsd_lig.xvg"), xvg("rmsf.xvg"), xvg("rg.xvg")
t_bb, bb = rb[:, 0], rb[:, 1] * 10  # nm -> Angstrom
t_lg, lg = rl[:, 0], rl[:, 1] * 10
t_rg, rgv = rg[:, 0] / 1000, rg[:, 1] * 10  # ps -> ns, nm -> A


def describe(name, t, y):
    h = len(y) // 2
    slope = np.polyfit(t[h:], y[h:], 1)[0]  # A per ns, second half only
    say(
        f"{name}: mean {y.mean():.2f} +/- {y.std():.2f} A | 2nd-half mean {y[h:].mean():.2f} A | "
        f"2nd-half drift {slope:+.3f} A/ns"
    )


describe("Protein backbone RMSD", t_bb, bb)
describe("Ligand RMSD          ", t_lg, lg)
describe("Radius of gyration   ", t_rg, rgv)
say("(drift near 0 means a plateau; a clearly positive drift means still relaxing)")

# ---------- 2. per-frame protein-ligand analysis ----------
u = mda.Universe(f"{d}/md.gro", f"{d}/md_fit.xtc")
other = u.select_atoms("not protein and not resname SOL HOH WAT TIP3 NA CL K")
names = sorted(set(other.resnames))
if not names:
    sys.exit("No ligand found in md.gro")
LIG = sys.argv[2] if len(sys.argv) > 2 else names[0]
say(f"\nLigand residue name used: {LIG}  (all non-water residues found: {names})")

lig = u.select_atoms(f"resname {LIG} and not name H*")
prot = u.select_atoms("protein and not name H*")
nfr = len(u.trajectory)
times, mind, ncont = [], [], []
res_count = {}
for ts in u.trajectory:
    pairs, dist = capped_distance(
        lig.positions, prot.positions, 6.0, box=ts.dimensions, return_distances=True
    )
    times.append(ts.time / 1000)
    mind.append(dist.min() if len(dist) else np.nan)
    close = pairs[dist <= 4.0]
    seen = set(prot.resindices[close[:, 1]])
    ncont.append(len(seen))
    for r in seen:
        res_count[r] = res_count.get(r, 0) + 1
times, mind, ncont = np.array(times), np.array(mind), np.array(ncont)

say(
    f"\nMin ligand-protein heavy-atom distance: mean {np.nanmean(mind):.2f} A, "
    f"max {np.nanmax(mind):.2f} A"
)
say(f"Protein residues within 4 A of ligand: mean {ncont.mean():.1f} per frame")
say("Top contact residues (fraction of frames within 4 A):")
for r, c in sorted(res_count.items(), key=lambda kv: -kv[1])[:12]:
    res = u.residues[r]
    say(f"  {res.resname}{res.resid:<4d} {c / nfr * 100:5.1f}%")

# ---------- 3. hydrogen bonds ----------
hb_time = None
try:
    h = HBA(
        u,
        between=["protein", f"resname {LIG}"],
        donors_sel="name N* O*",
        hydrogens_sel="name H*",
        acceptors_sel="name N* O*",
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=150,
    )
    h.run()
    hb = h.results.hbonds
    occ = {}
    for fr, di, hi, ai, *_ in hb:
        occ.setdefault((int(di), int(ai)), set()).add(int(fr))
    say("\nTop protein-ligand H-bonds (occupancy over all frames):")
    if not occ:
        say("  none found with 3.0 A / 150 deg criteria")
    for o, (di, ai) in sorted(
        ((len(v) / nfr * 100, k) for k, v in occ.items()), reverse=True
    )[:8]:
        A, B = u.atoms[di], u.atoms[ai]
        say(
            f"  {o:5.1f}%  {A.resname}{A.resid}:{A.name} -> {B.resname}{B.resid}:{B.name}"
        )
    cnt = np.zeros(nfr)
    for fr in hb[:, 0].astype(int):
        cnt[fr] += 1
    hb_time = cnt
except Exception:
    say("\nH-bond analysis FAILED (not 'no H-bonds'):")
    traceback.print_exc()

# ---------- 4. figure ----------
fig, ax = plt.subplots(2, 3, figsize=(15, 8))
ax[0, 0].plot(t_bb, bb, label="protein backbone")
ax[0, 0].plot(t_lg, lg, label="ligand")
ax[0, 0].set(xlabel="time (ns)", ylabel="RMSD (A)", title="RMSD")
ax[0, 0].legend()
ax[0, 1].plot(rf[:, 0], rf[:, 1] * 10)
ax[0, 1].set(xlabel="residue", ylabel="RMSF (A)", title="RMSF (C-alpha)")
ax[0, 2].plot(t_rg, rgv)
ax[0, 2].set(xlabel="time (ns)", ylabel="Rg (A)", title="Radius of gyration")
ax[1, 0].plot(times, mind)
ax[1, 0].set(
    xlabel="time (ns)", ylabel="distance (A)", title="Min ligand-protein distance"
)
ax[1, 1].plot(times, ncont)
ax[1, 1].set(
    xlabel="time (ns)", ylabel="residues", title="Residues within 4 A of ligand"
)
if hb_time is not None:
    ax[1, 2].plot(times, hb_time)
    ax[1, 2].set(xlabel="time (ns)", ylabel="count", title="Protein-ligand H-bonds")
else:
    ax[1, 2].axis("off")
plt.tight_layout()
plt.savefig(f"{d}/analysis.png", dpi=200)
open(f"{d}/summary.txt", "w").write("\n".join(out) + "\n")
print(f"\nSaved {d}/analysis.png and {d}/summary.txt")
