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
    filter_opt_data(collection, "02_curate-data/opt_records_to_remove.dat")


if __name__ == "__main__":
    # check_td()
    check_opt()
