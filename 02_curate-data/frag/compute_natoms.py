# compute and update database natoms column after manually altering the table

from tqdm import tqdm

from query import mol_from_smiles
from store import Store

s = Store("store.sqlite")

pairs = []
for m in tqdm(s.get_molecules(), total=s.get_sizehint()):
    mol = mol_from_smiles(m.smiles)
    pairs.append((mol.GetNumAtoms(), m.id))

s.cur.executemany("UPDATE molecules SET natoms = ? WHERE id = ?", pairs)
s.con.commit()
