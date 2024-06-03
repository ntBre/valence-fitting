import logging
from functools import partial
from multiprocessing import Pool

import click
from openff.toolkit import ForceField
from rdkit import Chem
from tqdm import tqdm

from store import DBForceField, DBMol, Match, Store, elements_to_bits
from utils import find_matches, mol_from_smiles

logger = logging.getLogger(__name__)

# fmt: off
PTABLE = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds ", "Rg ", "Cn ", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
# fmt: on


def into_params(ff) -> list[tuple[str, Chem.Mol]]:
    "Convert a ForceField into a sequence of pid:Mol pairs"
    return [
        (p.id, Chem.MolFromSmarts(p.smirks))
        for p in ff.get_parameter_handler("ProperTorsions").parameters
    ]


def load_want(filename):
    with open(filename) as inp:
        return {line.strip() for line in inp}


class Filter:
    def apply(self, mol: DBMol) -> bool:
        raise NotImplementedError()


class ElementFilter(Filter):
    def __init__(self, elements: int):
        self.mask = elements

    def apply(self, mol: DBMol) -> bool:
        return (mol.elements | self.mask) == self.mask


class NatomsFilter(Filter):
    def __init__(self, natoms: int):
        self.natoms = natoms

    def apply(self, mol: DBMol) -> bool:
        return mol.natoms <= self.natoms


def symbols_to_bits(symbols: list[str]) -> int:
    atomic_nums = [PTABLE.index(sym) for sym in symbols]
    return elements_to_bits(atomic_nums)


def parse_filters(filters: list[str]) -> list[Filter]:
    ret = list()
    for filt in filters:
        fields = filt.strip().split(":")
        match fields:
            case ["elements", arg]:
                atomic_symbols = [s.strip() for s in arg.split(",")]
                ret.append(ElementFilter(symbols_to_bits(atomic_symbols)))
            case ["natoms", natoms]:
                ret.append(NatomsFilter(int(natoms)))
            case _:
                raise NotImplementedError(f"unrecognized filter: `{fields}`")
    return ret


def inner(m: DBMol, params, filters) -> tuple[str, set[str]]:
    "Returns a SMILES and its matching parameter IDs"
    if not all((f.apply(m) for f in filters)):
        return "", set()
    mol = mol_from_smiles(m.smiles)
    res = set(find_matches(params, mol).values())
    if len(res) == 0:
        logger.warning(f"no matches found for {m.smiles}")
    return m.smiles, res


def _main(nprocs, chunk_size, filters, store_name, ffname, want, limit):
    s = Store(store_name)
    ff = ForceField(ffname)
    params = into_params(ff)

    pid_to_smirks = {
        p.id: p.smirks
        for p in ff.get_parameter_handler("ProperTorsions").parameters
    }

    res = dict()
    all_mols = [s for s in s.get_molecules(limit)]
    unmatched = 0
    with Pool(processes=nprocs) as p:
        for smiles, matches in tqdm(
            p.imap(
                partial(inner, params=params, filters=filters),
                all_mols,
                chunksize=chunk_size,
            ),
            total=len(all_mols),
        ):
            for pid in matches & want:
                if pid not in res:
                    res[pid] = Match(pid_to_smirks[pid], pid, list())
                res[pid].molecules.append(smiles)
            else:
                unmatched += 1

    logger.warning(f"{unmatched} SMILES not matching desired parameters")

    ret = list(res.values())
    s.insert_forcefield(DBForceField(ffname, ret))

    if logger.level > logging.WARNING:
        for mat in ret:
            print(mat.pid, len(mat.molecules))

    return s, ret


@click.command()
@click.option("--nprocs", "-n", default=8)
@click.option("--chunk-size", "-c", default=32)
@click.option("--filter", "-x", "filters", multiple=True)
@click.option("--store-name", "-s", default="store.sqlite")
@click.option("--ffname", "-f", default="openff-2.1.0.offxml")
@click.option("--target-params", "-t", default="want.params")
@click.option("--limit", "-l", hidden=True, default=None)
def main(
    nprocs, chunk_size, filters, store_name, ffname, target_params, limit
):
    filters = parse_filters(filters)
    want = load_want(target_params)
    _main(nprocs, chunk_size, filters, store_name, ffname, want, limit)


if __name__ == "__main__":
    main()
