import time
from typing import Iterator

import numpy as np
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import BRICS, Recap
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


def summary(mols: list[RDMol]) -> tuple[int, float, int]:
    "Compute statistics on streaming mols, printing at the end"
    num_atoms = []
    for mol in mols:
        num_atoms.append(mol.GetNumAtoms())
    mn = np.min(num_atoms)
    mean = np.mean(num_atoms)
    mx = np.max(num_atoms)
    return mn, mean, mx


def to_rdkit(mols: Iterator[Molecule]) -> Iterator[RDMol]:
    return (mol.to_rdkit() for mol in mols)


def recap(mols, minFragmentSize):
    ret = {}
    for mol in mols:
        leaves = Recap.RecapDecompose(
            mol.to_rdkit(), minFragmentSize=minFragmentSize
        ).GetLeaves()
        for smi, node in leaves.items():
            ret[smi] = node.mol
    return list(ret.values())


def brics(mols, minFragmentSize):
    # brics returns only molecules or only smiles, unlike recap for some
    # unknown reason, so reconstruct the smiles, mol pairs like recap to
    # deduplicate by smiles
    ret = {}
    for mol in mols:
        d = BRICS.BRICSDecompose(
            mol.to_rdkit(), minFragmentSize=minFragmentSize, returnMols=True
        )
        for mol in d:
            ret[Chem.MolToSmiles(mol)] = mol
    return list(ret.values())


mols = load_smiles("100.smi")

# draw_molecules("100.html", to_rdkit(mols))

print("Algo Mols Frags Min Mean Max Time")

for min_frag in [0, 2, 4, 8]:
    start = time.time()
    frags = recap(mols, min_frag)
    stop = time.time()
    draw_molecules(f"recap.{min_frag}.html", frags)
    mn, mean, mx = summary(frags)
    t = stop - start
    print(
        f"RECAP-{min_frag} {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}"
    )

for min_frag in [0, 2, 4, 8]:
    start = time.time()
    frags = brics(mols, min_frag)
    stop = time.time()
    draw_molecules(f"brics.{min_frag}.html", frags)
    mn, mean, mx = summary(frags)
    t = stop - start
    print(
        f"BRICS-{min_frag} {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}"
    )
