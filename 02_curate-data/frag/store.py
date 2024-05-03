import sqlite3
from itertools import chain
from multiprocessing import Pool

import click
from openff.toolkit import Molecule
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from rdkit import Chem
from tqdm import tqdm

from cut_compound import Compound


class Store:
    def __init__(self, filename="store.sqlite", nprocs=8):
        self.con = sqlite3.connect(filename)
        self.cur = self.con.cursor()
        self.cur.execute(
            """CREATE TABLE IF NOT EXISTS molecules (
            id integer primary key,
            smiles text unique
            )
            """
        )
        self.nprocs = nprocs

    def insert_molecule(self, smiles: str):
        "Insert a single SMILES into the database"
        self.cur.execute(
            "INSERT OR IGNORE INTO molecules (smiles) VALUES (?1)",
            (smiles,),
        )
        self.con.commit()

    def insert_molecules(self, smiles: list[str]):
        "Insert multiple SMILES into the database"
        self.cur.executemany(
            "INSERT OR IGNORE INTO molecules (smiles) VALUES (?1)",
            [(s,) for s in smiles],
        )
        self.con.commit()

    def process_line(line):
        [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
        try:
            mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
        except RadicalsNotSupportedError:
            return
        return xff(mol)

    def load_chembl(self, filename) -> dict[str, Molecule]:
        with open(filename) as inp, Pool(processes=self.nprocs) as pool:
            for frags in tqdm(
                pool.imap_unordered(
                    Store.process_line,
                    (line for i, line in enumerate(inp) if i > 0),
                    chunksize=32,
                ),
                total=2372675,
            ):
                if frags:
                    self.insert_molecules(frags)


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


def xff(mol):
    try:
        c = Compound(mol.to_rdkit())
    except Exception as e:
        print(f"warning: failed to convert to rdmol with {e}")
        return
    frags = c.cutCompound()

    ret = set()
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
                ret.add(Chem.MolToSmiles(frag))
        else:
            ret.add(Chem.MolToSmiles(emol))

    return ret


@click.command()
@click.option("--nprocs", "-n", type=int, default=8)
def main(nprocs):
    store = Store(nprocs=nprocs)
    store.load_chembl("chembl_33_chemreps.txt")


if __name__ == "__main__":
    main()
