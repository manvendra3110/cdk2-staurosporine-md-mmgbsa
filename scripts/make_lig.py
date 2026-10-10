import sys
from rdkit import Chem
from rdkit.Chem import AllChem

smi, pose, out = sys.argv[1:4]
m = Chem.MolFromPDBFile(pose, removeHs=True)
m = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smi), m)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, out)
print("net charge:", Chem.GetFormalCharge(m))
