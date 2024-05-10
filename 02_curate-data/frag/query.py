import json
from multiprocessing import Pool

import numpy as np
from openff.toolkit import ForceField
from rdkit import Chem
from tqdm import tqdm

from store import Store


def mol_from_smiles(smiles):
    """Create an RDKit molecule from SMILES and perform the cleaning operations
    from the OpenFF toolkit

    """
    rdmol = Chem.MolFromSmiles(smiles)
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


def find_smarts_matches(mol, smirks: Chem.Mol):
    "Adapted from RDKitToolkitWrapper._find_smarts_matches"
    idx_map = dict()
    for atom in smirks.GetAtoms():
        smirks_index = atom.GetAtomMapNum()
        if smirks_index != 0:
            idx_map[smirks_index - 1] = atom.GetIdx()
    map_list = [idx_map[x] for x in sorted(idx_map)]

    max_matches = np.iinfo(np.uintc).max
    match_kwargs = dict(
        uniquify=False, maxMatches=max_matches, useChirality=True
    )
    full_matches = mol.GetSubstructMatches(smirks, **match_kwargs)

    matches = [tuple(match[x] for x in map_list) for match in full_matches]

    return matches


def find_matches(
    params: list[tuple[str, Chem.Mol]], mol: Chem.Mol
) -> dict[tuple[int], str]:
    """Returns a map of chemical environment tuples to their matching parameter
    ids. Modeled (loosely) after ForceField.label_molecules and
    ParameterHandler._find_matches"""
    matches = {}
    for id, smirks in params:
        env_matches = find_smarts_matches(mol, smirks)
        for mat in env_matches:
            if mat[0] > mat[-1]:
                mat = mat[::-1]
            matches[mat] = id
    return matches


def into_params(ff) -> list[tuple[str, Chem.Mol]]:
    "Convert a ForceField into a sequence of pid:Mol pairs"
    return [
        (p.id, Chem.MolFromSmarts(p.smirks))
        for p in ff.get_parameter_handler("ProperTorsions").parameters
    ]


def load_want(filename):
    with open(filename) as inp:
        return {line.strip() for line in inp}


class Match:
    smirks: str
    pid: str
    molecules: list[str]

    def __init__(self, smirks, pid, molecules):
        self.smirks = smirks
        self.pid = pid
        self.molecules = molecules

    def to_dict(self):
        return dict(smirks=self.smirks, pid=self.pid, molecules=self.molecules)


def inner(smiles):
    mol = mol_from_smiles(smiles)
    return smiles, set(find_matches(params, mol).values())


if __name__ == "__main__":
    s = Store("store.sqlite")
    ff = ForceField(
        "../../01_generate-forcefield/output/initial-force-field-openff-2.1.0"
        ".offxml"
    )
    pid_to_smirks = {
        p.id: p.smirks
        for p in ff.get_parameter_handler("ProperTorsions").parameters
    }
    params = into_params(ff)
    want = load_want("want.params")
    res = dict()
    NPROCS = 8
    CHUNKSIZE = 32
    all_smiles = [s for s in s.get_smiles()]
    with Pool(processes=NPROCS) as p:
        for smiles, matches in tqdm(
            p.imap(inner, all_smiles, chunksize=CHUNKSIZE),
            total=len(all_smiles),
        ):
            for pid in matches & want:
                if pid not in res:
                    res[pid] = Match(pid_to_smirks[pid], pid, list())
                res[pid].molecules.append(smiles)

    with open("out.query.json", "w") as out:
        json.dump(res, out, default=Match.to_dict)

    # just looping smiles took ~7 seconds
    # MolFromSmiles took 4007 seconds (1:06:48)
    # mol_from_smiles took 5179 seconds (1:26:20)
