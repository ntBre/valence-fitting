import cProfile
from pstats import SortKey

from openff.toolkit import Molecule

from store import xff

with open("100.smi") as inp:
    smiles = [line.strip() for line in inp]

with cProfile.Profile() as pr:
    res = {
        f
        for mol in smiles
        for f in xff(Molecule.from_smiles(mol, allow_undefined_stereo=True))
    }
    pr.print_stats(SortKey.CUMULATIVE)

got = len(res)
want = 1009

assert got == want, f"{got} != {want}"
