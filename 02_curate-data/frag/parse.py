import json

from rdkit.Chem import AllChem

from query import Match, mol_from_smiles

with open("out.query.json") as f:
    data: dict[str, dict] = json.load(f)

matches: dict[str, Match] = {k: Match(**v) for k, v in data.items()}

fpgen = AllChem.GetMorganGenerator(radius=3)

for m in matches.values():
    lm = len(m.molecules)
    print(m.pid, m.smirks, lm)
    if lm <= 5:
        for smile in m.molecules:
            mol = mol_from_smiles(smile)
            print("\t", smile, len(fpgen.GetFingerprint(mol)))
