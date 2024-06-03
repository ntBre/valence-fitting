from openff.toolkit import ForceField

from query import into_params
from serve import ffname
from utils import find_matches, mol_from_smiles

ff = ForceField(ffname)
params = into_params(ff)
with open("want.smiles") as inp:
    for line in inp:
        if line.startswith("#"):
            continue
        pid, smirks, *smiles = line.split()
        for smile in smiles:
            mol = mol_from_smiles(smile)
            matches = find_matches(params, mol)
            if pid not in matches.values():
                for k, v in matches.items():
                    print(k, v)
                assert False, f"failed for {pid} ({smirks}) with {smile}"

print("all checks passed")
