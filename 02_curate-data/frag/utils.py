from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

from query import find_matches


def mol_from_smiles(smiles: str) -> Chem.Mol:
    """Create an RDKit molecule from SMILES and perform the cleaning operations
    from the OpenFF toolkit

    """
    rdmol = Chem.MolFromSmiles(smiles)
    return openff_clean(rdmol)


def mol_from_smarts(smarts: str) -> Chem.Mol:
    """Create an RDKit molecule from SMARTS and perform the cleaning operations
    from the OpenFF toolkit

    """
    rdmol = Chem.MolFromSmarts(smarts)
    return openff_clean(rdmol)


def openff_clean(rdmol: Chem.Mol) -> Chem.Mol:
    Chem.SanitizeMol(
        rdmol,
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
        ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    )
    Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
    Chem.AssignStereochemistry(rdmol)
    rdmol = Chem.AddHs(rdmol)
    return rdmol


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
