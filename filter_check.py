# checking to see which of our filters is removing all of the new torsion
# multiplicity data

import importlib
import logging
import sys
import typing

import numpy as np
from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    ElementFilter,
    HydrogenBondFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
)
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
from openff.toolkit.utils.toolkit_registry import _toolkit_registry_manager
from qcportal.models.records import RecordStatusEnum

sys.path.append("02_curate-data")
filters = importlib.import_module("filters")

logging.getLogger("openff").setLevel(logging.ERROR)


# copy pasta from 02/filters.py
def filter_td_data(
    dataset: "TorsionDriveResultCollection",
    td_records_to_remove: typing.Optional[str] = None,
    include_iodine: bool = False,
):
    "Filter a TorsionDrive dataset"

    if td_records_to_remove is not None:
        records_to_remove = np.loadtxt(td_records_to_remove, dtype=str)
    else:
        records_to_remove = []

    key = list(dataset.entries.keys())[0]

    print(f"initial len: {len(dataset.entries[key])}")
    # filter out entries to remove
    dataset.entries[key] = [
        entry
        for entry in dataset.entries[key]
        if entry.record_id not in records_to_remove
    ]
    print(f"after records to remove: {len(dataset.entries[key])}")

    # in a number of datasets the iodine-containing molecules
    # were tainted due to an auxiliary basis set issue
    # This has since been resolved and entries have been recomputed
    # in new datasets, but we still need to filter the old ones
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
    if include_iodine:
        elements.append("I")

    # filter out other unsuitable entries
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
        print(f"len {i}: {len(dataset.entries[key])}")

    return dataset


if __name__ == "__main__":
    with _toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper()])):
        collection = TorsionDriveResultCollection.parse_file(
            "02_curate-data/datasets/core-td.json"
        )
        filter_td_data(
            collection,
            "02_curate-data/td_records_to_remove.dat",
            False,
        )
