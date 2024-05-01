import sqlite3
from itertools import chain
from multiprocessing import Process, Queue

from openff.toolkit import Molecule
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from rdkit import Chem
from tqdm import tqdm

from cut_compound import Compound


class Store:
    def __init__(self, filename="store.sqlite"):
        self.con = sqlite3.connect(filename)
        self.cur = self.con.cursor()
        self.cur.execute(
            """CREATE TABLE IF NOT EXISTS molecules (
            id integer primary key,
            smiles text unique
            )
            """
        )
        # spawn a background process for inserting SMILES into the database
        # through a channel to avoid concurrent write attempts. calling
        # `self.send_molecules` enqueues a list of SMILES, which is dequeued in
        # the secondary process by `self.receive_molecules`, which calls the
        # synchronous `self.insert_molecules`
        self.queue = Queue()
        self.process = Process(target=self.receive_molecules)
        self.process.start()

    def insert_molecule(self, smiles: str):
        self.cur.execute(
            "INSERT OR IGNORE INTO molecules (smiles) VALUES (?1)",
            (smiles,),
        )
        self.con.commit()

    def insert_molecules(self, smiles: list[str]):
        self.cur.executemany(
            "INSERT OR IGNORE INTO molecules (smiles) VALUES (?1)",
            [(s,) for s in smiles],
        )
        self.con.commit()

    def receive_molecules(self):
        while True:
            self.insert_molecules(self.queue.get())

    def send_molecules(self, smiles: list[str]):
        self.queue.put(smiles)

    def load_chembl(self, filename, max_mols=1000) -> dict[str, Molecule]:
        found = 0
        with open(filename) as inp:
            for i, line in tqdm(enumerate(inp), total=2372675):
                if i == 0:  # skip header
                    continue
                [_chembl_id, cmiles, _inchi, _inchikey] = line.split("\t")
                try:
                    mol = Molecule.from_smiles(
                        cmiles, allow_undefined_stereo=True
                    )
                except RadicalsNotSupportedError:
                    continue
                frags = xff(mol)
                self.send_molecules(frags)
                found += len(frags)
                if max_mols and found >= max_mols:
                    break


def xff(mol):
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

    rdmol = mol.to_rdkit()
    c = Compound(rdmol)
    frags = c.cutCompound()

    ret = {}
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
                ret[Chem.MolToSmiles(frag)] = frag
        else:
            ret[Chem.MolToSmiles(emol)] = emol

    return ret


if __name__ == "__main__":
    store = Store()
    store.load_chembl("chembl_33_chemreps.txt", max_mols=None)
