import json

import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem

from query import Match, mol_from_smiles


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

    return ret


with open("out.query.json") as f:
    data: dict[str, dict] = json.load(f)

matches: dict[str, Match] = {k: Match(**v) for k, v in data.items()}

fpgen = AllChem.GetMorganGenerator(radius=3)

for m in matches.values():
    print(m.pid, m.smirks, len(m.molecules))
    fps = []
    for smile in m.molecules:
        mol = mol_from_smiles(smile)
        fps.append(fpgen.GetFingerprint(mol))
    dist = tanimoto(fps)
    print(dist)
