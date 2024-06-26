import re
from collections import defaultdict
from itertools import chain
from typing import Iterator

from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit import ForceField, Molecule
from tqdm import tqdm


def molecules(smiles: Iterator[str]) -> Iterator[Molecule]:
    seen = set()
    for s in smiles:
        if s not in seen:
            seen.add(s)
            yield Molecule.from_mapped_smiles(s, allow_undefined_stereo=True)


def bonds(mol: Molecule, i: int, j: int) -> int:
    "Return the multiplicity of the central bond between atoms i and j"

    def inner(mol, idx):
        return sum(
            (
                1
                for b in mol.bonds
                if b.atom1_index == idx or b.atom2_index == idx
            )
        )

    # -1 to discard shared central bond
    a = inner(mol, i) - 1
    b = inner(mol, j) - 1

    return a * b


def load_dataset(ds) -> Iterator[Molecule]:
    return molecules((e.cmiles for e in next(iter(ds.entries.values()))))


def param_sort_key(pid: str):
    m = re.match(r"t([0-9]+)(a-z)*", pid)
    return int(m[1]), m[2]


def main():
    # for original splits
    ff = ForceField("../openff-2.2.0-rc1.offxml")

    # for checking that splits worked - should output nothing
    # ff = ForceField("../output/initial-force-field-openff-2.1.0.offxml")

    opt = load_dataset(
        OptimizationResultCollection.parse_file(
            "../../02_curate-data/datasets/combined-opt.json"
        )
    )
    td = load_dataset(
        TorsionDriveResultCollection.parse_file(
            "../../02_curate-data/datasets/combined-td.json"
        )
    )
    bench = load_dataset(
        OptimizationResultCollection.parse_file(
            "../../../benchmarking/datasets/industry.json"
        )
    )
    mols = chain(opt, td, bench)

    def def_mol_dict():
        return defaultdict(set)

    res = defaultdict(set)
    mol_res = defaultdict(def_mol_dict)
    for mol in tqdm(mols, total=12475):
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        for (_, i, j, _), p in labels.items():
            mult = bonds(mol, i, j)
            res[p.id].add(mult)
            mol_res[p.id][mult].add(mol.to_smiles())

    h = ff.get_parameter_handler("ProperTorsions")
    for k, v in sorted(res.items(), key=lambda x: param_sort_key(x[0])):
        if len(v) > 1 and k not in [
            "t165",
            "t166",
            "t167",
        ]:  # skip k=0 torsions
            p = h.get_parameter(dict(id=k))[0]
            print(k, v, p.smirks)

            # v is a set of multiplicities, so loop over it and print up to the
            # 5 shortest smiles for each one. mol_res[k][i] is a set of smiles
            # corresponding to k (a parameter ID) and i (a single multiplicity)
            for i in v:
                smiles = list(sorted(mol_res[k][i], key=lambda s: len(s)))
                for smile in smiles[:5]:
                    print("\t", i, smile)


if __name__ == "__main__":
    main()
