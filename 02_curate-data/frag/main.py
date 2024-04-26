import time
from typing import Iterator

import matplotlib.pyplot as plt
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


class Summary:
    def __init__(self, mols: list[RDMol]):
        "Compute statistics on streaming mols, printing at the end"
        self.num_atoms = []
        for mol in mols:
            self.num_atoms.append(mol.GetNumAtoms())

    def stats(self) -> tuple[int, float, int]:
        mn = np.min(self.num_atoms)
        mean = np.mean(self.num_atoms)
        mx = np.max(self.num_atoms)
        return mn, mean, mx

    def hist(self, ax, nbins, title):
        ax.hist(
            self.num_atoms, nbins, histtype="step", stacked=True, fill=False
        )
        ax.set_title(title)


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


def erb(mols):
    "Fragment on every rotatable bond, adapted from Lily's script"
    ret = {}
    for mol in mols:
        bonds = mol.find_rotatable_bonds()
        rdmol = mol.to_rdkit()
        rdbonds = tuple(
            [
                rdmol.GetBondBetweenAtoms(
                    bond.atom1_index, bond.atom2_index
                ).GetIdx()
                for bond in bonds
            ]
        )
        fragmented = Chem.FragmentOnBonds(rdmol, rdbonds)
        for frag in Chem.GetMolFrags(fragmented, asMols=True):
            ret[Chem.MolToSmiles(frag)] = frag
    return list(ret.values())


mols = load_smiles("100.smi")

# draw_molecules("100.html", to_rdkit(mols))

fig, ((ax0, ax1, ax2), (ax3, ax4, ax5), (ax6, ax7, ax8)) = plt.subplots(3, 3)
axes = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
ptr = 0

NBINS = 10

print("Algo Mols Frags Min Mean Max Time")

for mf in [0, 2, 4, 8]:
    start = time.time()
    frags = recap(mols, mf)
    stop = time.time()
    draw_molecules(f"output/recap.{mf}.html", frags)
    s = Summary(frags)
    mn, mean, mx = s.stats()
    title = f"RECAP-{mf}"
    s.hist(axes[ptr], NBINS, title)
    ptr += 1
    t = stop - start
    print(f"{title} {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}")

for mf in [0, 2, 4, 8]:
    start = time.time()
    frags = brics(mols, mf)
    stop = time.time()
    draw_molecules(f"output/brics.{mf}.html", frags)
    s = Summary(frags)
    mn, mean, mx = s.stats()
    title = f"BRICS-{mf}"
    s.hist(axes[ptr], NBINS, title)
    ptr += 1
    t = stop - start
    print(f"{title} {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}")

start = time.time()
frags = erb(mols)
stop = time.time()
draw_molecules("output/erb.html", frags)
s = Summary(frags)
mn, mean, mx = s.stats()
title = "ERB"
s.hist(axes[ptr], NBINS, title)
t = stop - start
print(f"ERB {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}")

fig.tight_layout()
plt.savefig("hist.png")
