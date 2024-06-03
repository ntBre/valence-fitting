from openff.toolkit import ForceField

from query import find_matches, into_params, load_want, mol_from_smiles
from store import Store

ffname = (
    "../../01_generate-forcefield/output/"
    "initial-force-field-openff-2.1.0.offxml"
)
ff = ForceField(ffname)
params = into_params(ff)
s = Store()

want = load_want("want.params")

print(want)

for m in s.get_molecules(limit=500):
    mol = mol_from_smiles(m.smiles)
    res = set(find_matches(params, mol).values())
    if m.natoms < 25:
        print(m.smiles, m.natoms, mol.GetNumAtoms(), res)
