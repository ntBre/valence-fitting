from typing import Iterator

import numpy as np
from openff.toolkit import Molecule
from rdkit.Chem import Recap
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D
from rdkit.Chem.rdchem import Mol as RDMol


def load_chembl(filename, max_atoms=80, max_mols=1000) -> list[Molecule]:
    ret = []
    with open(filename) as inp:
        for i, line in enumerate(inp):
            if i == 0:  # skip header
                continue
            [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
            mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
            if mol.n_atoms < max_atoms:
                ret.append(mol)
            if len(ret) >= max_mols:
                break
    return ret


def save_mols(filename, mols):
    with open(filename, "w") as out:
        for mol in mols:
            print(mol.to_smiles(), file=out)


def load_smiles(filename) -> list[Molecule]:
    with open(filename) as inp:
        return [
            Molecule.from_smiles(s.strip(), allow_undefined_stereo=True)
            for s in inp
        ]


def draw_rdkit(rdmol) -> str:
    "Return a single molecule drawn as an SVG string"
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(rdmol)
    rdmol = rdMolDraw2D.PrepareMolForDrawing(rdmol)
    return MolsToGridImage(
        [rdmol],
        useSVG=True,
        subImgSize=(300, 300),
        molsPerRow=1,
    )


def draw_molecules(filename, mols: Iterator[RDMol]):
    with open(filename, "w") as out:
        for mol in mols:
            print(draw_rdkit(mol), file=out)


def summary(mols: Iterator[RDMol]) -> Iterator[RDMol]:
    "Compute statistics on streaming mols, printing at the end"
    num_atoms = []
    for mol in mols:
        num_atoms.append(mol.GetNumAtoms())
        yield mol
    print(f"{len(num_atoms)} molecules")
    mn = np.min(num_atoms)
    mean = np.mean(num_atoms)
    mx = np.max(num_atoms)
    print(f"{mn}, {mean:.2f}, {mx} (min, mean, max) atoms")


def to_rdkit(mols: Iterator[Molecule]) -> Iterator[RDMol]:
    return (mol.to_rdkit() for mol in mols)


def recap(mols, minFragmentSize):
    ret = {}
    for mol in mols:
        rdmol = mol.to_rdkit()
        leaves = Recap.RecapDecompose(
            rdmol, minFragmentSize=minFragmentSize
        ).GetLeaves()
        for smi, node in leaves.items():
            ret[smi] = node.mol
    return list(ret.values())


mols = load_smiles("100.smi")

# draw_molecules("100.html", to_rdkit(mols))

for min_frag in [0, 2, 4, 8]:
    print(f"Recap {min_frag}")
    draw_molecules(f"recap.{min_frag}.html", summary(recap(mols, min_frag)))
    print("=======")
