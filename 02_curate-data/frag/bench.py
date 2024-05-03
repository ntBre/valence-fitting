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

# Conclusion: all of the time is spent in xff as expected, and most of the time
# in that is spent computing graph isomorphisms in networkx, so it's not clear
# to me how to speed anything up. my minor changes had no effect on the
# runtime, but hopefully removing the global variable avoids any potential
# issues with using multithreading. on the other hand, setting the chain_num to
# zero in every case also didn't change the final result length here either.
