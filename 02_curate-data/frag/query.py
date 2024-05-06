import time

from rdkit import Chem
from tqdm import tqdm

from store import Store


def mol_from_smiles(smiles):
    """Create an RDKit molecule from SMILES and perform the cleaning operations
    from the OpenFF toolkit

    """
    rdmol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(
        rdmol,
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    )
    Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
    Chem.AssignStereochemistry(rdmol)
    Chem.AddHs(rdmol)
    return rdmol


s = Store("store.sqlite")

start = time.time()

i = 0
for smiles in tqdm(s.get_smiles(), total=s.get_sizehint()):
    mol_from_smiles(smiles)
    i += 1

end = time.time()

# just looping smiles took ~7 seconds
# MolFromSmiles took 4007 seconds (1:06:48)
# mol_from_smiles took 5179 seconds (1:26:20)
print(f"processed {i} molecules in {end-start:.2f} sec")
