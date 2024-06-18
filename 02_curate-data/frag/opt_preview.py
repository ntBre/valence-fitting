import logging
import re

from openff.toolkit import ForceField, Molecule
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

from serve import ffname

logging.getLogger("openff").setLevel(logging.ERROR)


def mol_to_svg(mol: Molecule, hl_atoms) -> str:
    mol = mol.to_rdkit()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    rdmol = rdMolDraw2D.PrepareMolForDrawing(mol)

    return MolsToGridImage(
        [rdmol],
        useSVG=True,
        highlightAtomLists=[hl_atoms],
        subImgSize=(300, 300),
        molsPerRow=1,
    )


junk = re.compile("[(,)]")

ff = ForceField(ffname)

with open("opt.smi") as inp, open("opt.html", "w") as out:
    for line in inp:
        if line.startswith("#"):
            continue
        pid, cmiles, *_ = line.split()
        mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        d = {p.id: env for env, p in labels.items()}
        assert pid in d

        print(
            '<span style="display:inline-block;border:1px solid black">',
            file=out,
        )
        print(f"<h1>{pid}</h1>", file=out)
        print(mol_to_svg(mol, d[pid]), file=out)
        print("</span>", file=out)
