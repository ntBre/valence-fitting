import time

from rdkit import Chem
from tqdm import tqdm

from store import Store

s = Store("store.sqlite")

start = time.time()

i = 0
for smiles in tqdm(s.get_smiles(), total=s.get_sizehint()):
    Chem.MolFromSmiles(smiles)
    i += 1

end = time.time()

# just looping smiles took ~7 seconds
print(f"processed {i} molecules in {end-start:.2f} sec")
