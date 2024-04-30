import time
import warnings
from itertools import chain
from typing import Iterator

import matplotlib.pyplot as plt
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import BRICS, Recap
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D
from rdkit.Chem.rdchem import Mol as RDMol

from cut_compound import Compound

warnings.simplefilter("ignore")
with warnings.catch_warnings():
    import numpy as np


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
    print(Chem.MolToSmiles(rdmol))
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


def xff(mols):
    def find_frag_bonds(rdmol, keep_atoms):
        "Locate bonds between atoms to keep and those to remove"
        for bond in rdmol.GetBonds():
            print(bond)
        panic
    ret = {}
    for mol in mols:
        rdmol = mol.to_rdkit()
        c = Compound(rdmol)
        frags = c.cutCompound()
        all_atoms = set(range(c.rdmol.GetNumAtoms()))
        for fr in chain(frags.frag_rings, frags.frag_chains):
            emol = Chem.rdchem.EditableMol(c.rdmol)
            frag_atoms = set(fr)
            to_remove = all_atoms - frag_atoms
            for idx in sorted(to_remove, reverse=True):
                emol.RemoveAtom(idx)
            mol = emol.GetMol()
            ret[Chem.MolToSmiles(mol)] = mol
    return list(ret.values())


def run_algo(fun, mols, html, title, *args):
    global ptr
    start = time.time()
    frags = fun(mols, *args)
    stop = time.time()
    draw_molecules(html, frags)
    s = Summary(frags)
    mn, mean, mx = s.stats()
    s.hist(axes[ptr], NBINS, title)
    ptr += 1
    t = stop - start
    print(f"{title} {len(mols)} {len(frags)} {mn} {mean:.2f} {mx} {t:.1f}")


mols = load_smiles("100.smi")

# draw_molecules("100.html", to_rdkit(mols))

fig, ((ax0, ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8, ax9)) = plt.subplots(
    2, 5
)
axes = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
ptr = 0

NBINS = 10

print("Algo Mols Frags Min Mean Max Time")

# for mf in [0, 2, 4, 8]:
#     run_algo(recap, mols, f"output/recap.{mf}.html", f"RECAP-{mf}", mf)

# for mf in [0, 2, 4, 8]:
#     run_algo(brics, mols, f"output/brics.{mf}.html", f"BRICS-{mf}", mf)

# run_algo(erb, mols, "output/erb.html", "ERB")

import faulthandler

with open("fault_handler.log", "w") as f:
    faulthandler.enable(f)

    run_algo(xff, mols, "output/xff.html", "XFF")

fig.tight_layout()
plt.savefig("hist.png")
