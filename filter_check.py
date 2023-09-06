# checking to see which of our filters is removing all of the new torsion
# multiplicity data

import importlib
import logging
import sys
import warnings

from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.qcsubmit.results.filters import (
    ConformerRMSDFilter,
    ConnectivityFilter,
    ElementFilter,
    HydrogenBondFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
)
from qcportal.models.records import RecordStatusEnum
from tqdm import tqdm

sys.path.append("02_curate-data")
filters = importlib.import_module("filters")

logging.getLogger("openff").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning)


# copy pasta from 02/filters.py
def filter_td_data(dataset: "TorsionDriveResultCollection"):
    key = list(dataset.entries.keys())[0]
    last = len(dataset.entries[key])
    print(f"initial len: {last}")
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
    for i, filt in enumerate(
        [
            HydrogenBondFilter(method="baker-hubbard"),
            RecordStatusFilter(status=RecordStatusEnum.complete),
            ConnectivityFilter(tolerance=1.2),
            UnperceivableStereoFilter(),
            ElementFilter(allowed_elements=elements),
            filters.ChargeCheckFilter(),
        ]
    ):
        dataset = dataset.filter(filt)
        new = len(dataset.entries[key])
        print(f"len after {filt.__class__}: {new} Δ: {last - new}")
        last = new
    return dataset


def filter_opt_data(dataset):
    key = list(dataset.entries.keys())[0]
    last = len(dataset.entries[key])
    print(f"initial len: {last}")
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
    for filt in [
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        ElementFilter(allowed_elements=elements),
        ConformerRMSDFilter(max_conformers=12),
        filters.ChargeCheckFilter(),
    ]:
        dataset = dataset.filter(filt)
        new = len(dataset.entries[key])
        print(f"len after {filt.__class__}: {new} Δ: {last - new}")
        last = new

    return dataset


def check_td():
    collection = TorsionDriveResultCollection.parse_file(
        "02_curate-data/datasets/core-td.json"
    )
    filter_td_data(collection)


def check_opt():
    collection = OptimizationResultCollection.parse_file(
        "02_curate-data/datasets/core-opt.json"
    )
    filter_opt_data(collection)


if __name__ == "__main__":
    # check_td()
    # check_opt()

    # I'm getting a mismatch between the number from check_opt and when I ran
    # the filter initially. the check is giving me 25 from connectivity + 119
    # from conformers + 38 from the charge check = 182 filtered out = 218
    # remaining. in contrast, all of the ways I'm checking how many should be
    # left (grep record_id in filtered-opt.json, composition.py, fast-filter
    # output) all give 257 remaining or 143 filtered out. this gives 218 too.
    # is there something wrong with fast-filter??

    # there is indeed something wrong with fast-filter. passing a batch size of
    # 400 (to get everything in one python file) gives the expected pure python
    # result of 218. clearly something is going wrong with the batching

    # let's see if batching causes the same issue in Python proper

    # doing the batching in Python gives the same result!

    # what about using the filter_opt_data function in this file?

    # same result using batching here!

    # next question is whether a batch size of 400 restores 218 like it did for
    # fast-filter. that will rule out an issue in fast-filter completely. I'm
    # basically trying to assemble a minimal example to report

    # batch size 400 goes back to 218, time for a bug report

    collection = OptimizationResultCollection.parse_file(
        "02_curate-data/datasets/core-opt.json"
    )

    batch_size = 400

    key = list(collection.entries.keys())[0]
    entries = [v for vals in collection.entries.values() for v in vals]

    batches = [
        OptimizationResultCollection(
            entries={key: entries[i : min(i + batch_size, len(entries))]}
        )
        for i in range(0, len(entries), batch_size)
    ]

    final = []
    for batch in tqdm(batches):
        batch = filter_opt_data(batch)
        final.extend(v for vals in batch.entries.values() for v in vals)

    print(len(final))
