from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

from query import find_matches


def make_svg(pid, map, mol_map, mol):
    hl_atoms = []
    if pid in map:
        tmp = find_matches(mol_map, mol)
        hl_atoms, _pid = [
            (atoms, param_id)
            for atoms, param_id in tmp.items()
            if param_id == pid
        ][0]
    return mol_to_svg(mol, 400, 300, "", hl_atoms)


def mol_to_svg(
    mol: Chem.Mol, width: int, height: int, legend: str, hl_atoms: list[int]
) -> str:
    "Return an SVG of `mol`"
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


def find_smallest(mols: list[Chem.Mol], cluster: list[int]) -> int:
    "Return the index of the smallest molecule in `cluster`"
    tmp = [(i, mols[c]) for i, c in enumerate(cluster)]
    return min(tmp, key=lambda x: x[1].GetNumAtoms())[0]
