from itertools import chain

from openff.toolkit import Molecule
from rdkit import Chem
from tqdm import tqdm

from cut_compound import Compound


def xff(mol):
    def find_frag_bonds(rdmol, keep_atoms):
        "Locate bonds between atoms to keep and those to remove"
        to_remove = []
        for bond in rdmol.GetBonds():
            b1, b2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if (b1 in keep_atoms and b2 not in keep_atoms) or (
                b2 in keep_atoms and b1 not in keep_atoms
            ):
                to_remove.append(rdmol.GetBondBetweenAtoms(b1, b2).GetIdx())
        return to_remove

    rdmol = mol.to_rdkit()
    c = Compound(rdmol)
    frags = c.cutCompound()

    ret = {}
    for fr in chain(frags.frag_rings, frags.frag_chains):
        emol = Chem.Mol(c.rdmol)
        frag_atoms = set(fr)
        rdbonds = find_frag_bonds(emol, frag_atoms)
        if rdbonds:
            fragmented = Chem.FragmentOnBonds(
                emol,
                rdbonds,
                addDummies=True,
                dummyLabels=[(0, 0)] * len(rdbonds),
            )
            for frag in Chem.GetMolFrags(fragmented, asMols=True):
                ret[Chem.MolToSmiles(frag)] = frag
        else:
            ret[Chem.MolToSmiles(emol)] = emol

    return ret


def load_chembl(filename, max_mols=1000) -> dict[str, Molecule]:
    ret = {}
    with open(filename) as inp:
        for i, line in tqdm(enumerate(inp)):
            if i == 0:  # skip header
                continue
            [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
            mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
            ret.update(xff(mol))
            if len(ret) >= max_mols:
                break
    return ret


if __name__ == "__main__":
    mols = load_chembl("chembl_33_chemreps.txt")
    for smiles in mols.keys():
        print(smiles)
