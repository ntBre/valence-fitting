from rdkit import Chem
from tqdm import tqdm

from store import Store
from utils import mol_from_smiles

s = Store()
res = []
for m in tqdm(s.get_molecules(), desc="Adding inchis", total=s.get_sizehint()):
    mol = mol_from_smiles(m.smiles)
    res.append((m.id, Chem.MolToInchiKey(mol)))

s.cur.executemany("UPDATE molecules SET inchikey = ?1 WHERE id = ?2", res)
s.con.commit()
