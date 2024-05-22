import sqlite3
from itertools import chain
from multiprocessing import Pool
from typing import Iterator

import click
from openff.toolkit import Molecule
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from rdkit import Chem
from tqdm import tqdm

from cut_compound import Compound


class DBMol:
    "Representation of a molecule in the database"
    id: int
    smiles: str
    natoms: int

    def __init__(self, id, smiles, natoms):
        self.id = id
        self.smiles = smiles
        self.natoms = natoms

    def __repr__(self):
        return (
            f"DBMol(id={self.id}, smiles={self.smiles}, natoms={self.natoms})"
        )


class Store:
    def __init__(self, filename="store.sqlite", nprocs=8):
        self.con = sqlite3.connect(filename)
        self.cur = self.con.cursor()
        self.cur.arraysize = 1024
        self.cur.execute(
            """CREATE TABLE IF NOT EXISTS molecules (
            id integer primary key,
            smiles text unique,
            natoms int
            )
            """
        )
        self.nprocs = nprocs

    def insert_molecule(self, smiles: str):
        "Insert a single SMILES into the database"
        self.insert_molecules([smiles])

    def insert_molecules(self, smiles: dict[str, int]):
        "Insert multiple SMILES into the database"
        self.cur.executemany(
            "INSERT OR IGNORE INTO molecules (smiles, natoms) VALUES (?1, ?2)",
            [(s, n) for s, n in smiles.items()],
        )
        self.con.commit()

    def get_sizehint(self) -> int:
        "Return a count of rows in the database"
        return self.cur.execute("SELECT COUNT(*) FROM molecules").fetchone()[0]

    def get_smiles(self, limit=None) -> Iterator[str]:
        "Return an iterator over SMILES in the database"
        if limit:
            res = self.cur.execute(
                "SELECT smiles FROM molecules limit ?", (limit,)
            )
        else:
            res = self.cur.execute("SELECT smiles FROM molecules")
        while len(v := res.fetchmany()) > 0:
            # unpack 1-tuples
            yield from (x[0] for x in v)

    def get_molecules(self) -> Iterator[DBMol]:
        res = self.cur.execute("SELECT id, smiles, natoms FROM molecules")
        while len(v := res.fetchmany()) > 0:
            # unpack 3-tuples
            yield from (DBMol(id=x[0], smiles=x[1], natoms=x[2]) for x in v)

    def process_line(line) -> dict[str, int]:
        [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
        all_smiles = cmiles.split(".")
        ret = dict()
        for smiles in all_smiles:
            try:
                mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
            except RadicalsNotSupportedError:
                continue
            if x := xff(mol):
                ret.update(x)
        return ret

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


def xff(mol) -> dict[str, int] | None:
    try:
        c = Compound(mol.to_rdkit())
    except Exception as e:
        print(f"warning: failed to convert to rdmol with {e}")
        return
    frags = c.cutCompound()

    ret = dict()
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
                ret[Chem.MolToSmiles(frag)] = frag.GetNumAtoms()
        else:
            ret[Chem.MolToSmiles(emol)] = emol.GetNumAtoms()

    return ret


@click.command()
@click.option("--nprocs", "-n", type=int, default=8)
def main(nprocs):
    store = Store(nprocs=nprocs)
    store.load_chembl("chembl_33_chemreps.txt")


if __name__ == "__main__":
    main()
