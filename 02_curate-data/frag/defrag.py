from openff.toolkit import ForceField
from rdkit import Chem

from query import load_want, mol_from_smarts
from serve import ffname
from store import Store

s = Store("store.sqlite")
want = load_want("want.params")
ff = ForceField(ffname)

smirks = {
    p.id: p.smirks
    for p in ff.get_parameter_handler("ProperTorsions").parameters
}

for pid in want:
    print(pid, smirks[pid])

# this doesn't work very well because the ~ bonds and stuff aren't filled in
# with actual bonds, and then AddHs overfills the valence of the atoms.
# Instead, I'm going to try to craft SMILES by hand that match each parameter
# with manual * insertions and then fill those with fragments from the database

for smirk in smirks.values():
    mol = Chem.MolFromSmarts(smirk)
    print(smirk, Chem.MolToSmiles(mol))
