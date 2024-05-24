import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem

from query import mol_from_smiles
from store import Store


def tanimoto(fps):
    "Compute the tanimoto distance matrix for a list of fingerprints"
    n = len(fps)
    ret = np.zeros((n, n))
    for row in range(n):
        ret[row, row:] = np.array(
            DataStructs.BulkTanimotoSimilarity(
                fps[row], fps[row:], returnDistance=True
            )
        )

    return ret + ret.T  # copy upper triangle to lower, diag is already zero


ff = (
    "../../01_generate-forcefield/output/"
    "initial-force-field-openff-2.1.0.offxml"
)
s = Store("store.sqlite")
dbff = s.get_forcefield(ff)

matches = dbff.matches

print(f"found {len(matches)} matches")

fpgen = AllChem.GetMorganGenerator(radius=3)

for m in matches:
    print(m.pid, m.smirks, len(m.molecules))
    fps = []
    for smile in m.molecules:
        mol = mol_from_smiles(smile)
        fps.append(fpgen.GetFingerprint(mol))
    dist = tanimoto(fps)
    print(dist)
