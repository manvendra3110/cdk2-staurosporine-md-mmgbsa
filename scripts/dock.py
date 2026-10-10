import csv, subprocess, os, re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from scipy.stats import spearmanr

xyz = np.array(
    [
        [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        for l in open("ligand_crystal.pdb")
        if l.startswith("HETATM")
    ]
)
cx, cy, cz = xyz.mean(axis=0)
os.makedirs("dock", exist_ok=True)


def load(pdb, smi):
    t = Chem.MolFromSmiles(smi)
    m = Chem.MolFromPDBFile(pdb, removeHs=True)
    return AllChem.AssignBondOrdersFromTemplate(t, m)


rows = list(csv.DictReader(open("ligands.csv")))
res = []
for r in rows:
    n, smi = r["name"], r["smiles"]
    q = f"dock/{n}.pdbqt"
    o = f"dock/{n}_out.pdbqt"
    subprocess.run(["obabel", f"-:{smi}", "-O", q, "--gen3d", "-p", "7.4"], check=True)
    subprocess.run(
        [
            "vina",
            "--receptor",
            "rec.pdbqt",
            "--ligand",
            q,
            "--center_x",
            str(cx),
            "--center_y",
            str(cy),
            "--center_z",
            str(cz),
            "--size_x",
            "25",
            "--size_y",
            "25",
            "--size_z",
            "25",
            "--exhaustiveness",
            "16",
            "--num_modes",
            "9",
            "--seed",
            "42",
            "--out",
            o,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    score = float(re.search(r"VINA RESULT:\s+(-?\d+\.\d+)", open(o).read()).group(1))
    subprocess.run(
        ["obabel", o, "-O", f"dock/{n}_pose.pdb", "-f", "1", "-l", "1"], check=True
    )
    rmsd = ""
    if n == "native":
        rmsd = rdMolAlign.CalcRMS(
            load(f"dock/{n}_pose.pdb", smi), load("ligand_crystal.pdb", smi)
        )
        print(f"REDOCKING RMSD = {rmsd:.2f} A  (goal: under 2.0)")
    res.append((n, r["ic50_nM"], score, rmsd))
    print(n, score)

with open("dock_results.csv", "w") as f:
    f.write("name,ic50_nM,vina_kcal,redock_rmsd\n")
    for x in res:
        f.write(",".join(map(str, x)) + "\n")

pairs = [(9 - np.log10(float(i)), -s) for n, i, s, _ in res if i]
if len(pairs) > 2:
    rho, p = spearmanr(*zip(*pairs))
    print(f"Docking score vs pIC50: Spearman rho = {rho:.2f}, p = {p:.3f}")
