# add the inchikeys from the new dataset (dataset.smi) to inchis.dat

from openff.toolkit import Molecule

with open("dataset.smi") as inp:
    for line in inp:
        if line.startswith("#"):
            continue
        pid, cmiles, *rest = line.split()
        mol = Molecule.from_mapped_smiles(cmiles, allow_undefined_stereo=True)
        print(mol.to_inchikey())
