from openff.toolkit import Molecule


def load_smiles(filename) -> list[Molecule]:
    ret = []
    with open(filename) as inp:
        for i, line in enumerate(inp):
            if i == 0:  # skip header
                continue
            [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
            mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
            if mol.n_atoms < 80:
                ret.append(mol)
            if len(ret) >= 1000:
                break
    return ret


def save_mols(filename, mols):
    with open(filename, "w") as out:
        for mol in mols:
            print(mol.to_smiles(), file=out)


save_mols("1000.smi", load_smiles("chembl_33_chemreps.txt"))
