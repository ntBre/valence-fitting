import json
import logging
import time
from collections import defaultdict
from typing import Union

import click
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit import ForceField, Molecule
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)


# copy pasta from known-issues
class Timer:
    def __init__(self):
        self._start = time.time()

    def say(self, s: str):
        now = time.time()
        print(f"{s} after {now - self._start:.1f} s")


# also copy pasta from known-issues
def load_dataset(
    dataset: str,
) -> Union[OptimizationResultCollection, TorsionDriveResultCollection]:
    """Peeks at the first entry of `dataset` to determine its type and
    then loads it appropriately.

    Raises a `TypeError` if the first entry is neither a `torsion`
    record nor an `optimization` record.
    """
    with open(dataset, "r") as f:
        j = json.load(f)
    entries = j["entries"]
    keys = entries.keys()
    assert len(keys) == 1  # only handling this case for now
    key = list(keys)[0]
    match j["entries"][key][0]["type"]:
        case "torsion":
            return TorsionDriveResultCollection.parse_file(dataset)
        case "optimization":
            return OptimizationResultCollection.parse_file(dataset)
        case t:
            raise TypeError(f"Unknown result collection type: {t}")


# helper class for returning a sized iterator
class Molecules:
    def __init__(self, dataset):
        self.data = [v for value in dataset.entries.values() for v in value]
        self.length = len(self.data)

    def __len__(self):
        return self.length

    def __iter__(self):
        for r in self.data:
            yield Molecule.from_mapped_smiles(
                r.cmiles, allow_undefined_stereo=True
            )


def to_molecules(dataset) -> Molecules:
    return Molecules(dataset)


# monkey patch for cute calls
TorsionDriveResultCollection.to_molecules = to_molecules
OptimizationResultCollection.to_molecules = to_molecules


@click.command()
@click.option("--forcefield", "-f", default="openff-2.1.0.offxml")
@click.option(
    "--dataset", "-d", default="02_curate-data/datasets/combined-td.json"
)
def check_coverage(forcefield, dataset):
    "Check proper torsion parameter coverage in `forcefield` using `dataset`"

    print("checking coverage with")
    print(f"forcefield = {forcefield}")
    print(f"dataset = {dataset}")

    timer = Timer()

    ff = ForceField(forcefield, allow_cosmetic_attributes=True)
    td_data = load_dataset(dataset)

    timer.say("finished loading collection")

    h = ff.get_parameter_handler("ProperTorsions")
    tors_ids = [p.id for p in h.parameters]

    records_and_molecules = td_data.to_molecules()

    timer.say("finished to_records")

    results = defaultdict(int)
    # max of 1 per molecule
    mol_results = defaultdict(int)
    for molecule in tqdm(records_and_molecules, desc="Counting results"):
        all_labels = ff.label_molecules(molecule.to_topology())[0]
        torsions = all_labels["ProperTorsions"]
        for torsion in torsions.values():
            results[torsion.id] += 1
        for tid in {t.id for t in torsions.values()}:
            mol_results[tid] += 1

    timer.say("finished counting results")

    got = len(results)
    want = len(tors_ids)
    pct = 100.0 * float(got) / float(want)
    print(f"{got} / {want} ({pct:.1f}%) tors ids covered:")

    for id in tors_ids:
        smirk = h.get_parameter(dict(id=id))[0].smirks
        print(f"{id:5}{results[id]:5}{mol_results[id]:5}   {smirk}")

    missing_ids = [k for k in results.keys() if results[k] == 0]
    missing_smirks = [
        h.get_parameter(dict(id=p))[0].smirks for p in missing_ids
    ]
    print("\nmissing ids:")
    for i, (id, smirk) in enumerate(zip(missing_ids, missing_smirks)):
        print(f"{i:5}{id:>7}   {smirk}")

    timer.say("finished")


if __name__ == "__main__":
    check_coverage()
