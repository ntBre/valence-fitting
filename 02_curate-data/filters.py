import logging
import typing

import numpy as np
from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConformerRMSDFilter,
    ConnectivityFilter,
    ElementFilter,
    HydrogenBondFilter,
    RecordStatusFilter,
    SinglepointRecordFilter,
    UnperceivableStereoFilter,
)
from openff.toolkit.utils.exceptions import (
    ChargeCalculationError,
    ConformerGenerationError,
)
from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper
from qcportal.record_models import RecordStatusEnum

logging.getLogger("openff").setLevel(logging.ERROR)


class ChargeCheckFilter(SinglepointRecordFilter):
    def _filter_function(self, result, record, molecule) -> bool:
        try:
            OpenEyeToolkitWrapper().assign_partial_charges(
                molecule, partial_charge_method="am1bccelf10"
            )
        except (ChargeCalculationError, ConformerGenerationError):
            return False
        else:
            return True


def filter_opt_data(
    dataset,
    opt_records_to_remove: typing.Optional[str] = None,
    include_iodine: bool = False,
    max_opt_conformers: int = 12,
):
    if opt_records_to_remove is not None:
        records_to_remove = np.loadtxt(opt_records_to_remove, dtype=str)
    else:
        records_to_remove = []

    key = list(dataset.entries.keys())[0]

    # filter out entries to remove
    dataset.entries[key] = [
        entry
        for entry in dataset.entries[key]
        if entry.record_id not in records_to_remove
    ]

    # in a number of datasets the iodine-containing molecules
    # were tainted due to an auxiliary basis set issue
    # This has since been resolved and entries have been recomputed
    # in new datasets, but we still need to filter the old ones
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
    if include_iodine:
        elements.append("I")

    # filter out other unsuitable entries
    dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        # UnperceivableStereoFilter(),
        # ElementFilter(allowed_elements=elements),
        # ConformerRMSDFilter(max_conformers=max_opt_conformers),
        ChargeCheckFilter(),
    )

    return dataset


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

    # filter out entries to remove
    dataset.entries[key] = [
        entry
        for entry in dataset.entries[key]
        if entry.record_id not in records_to_remove
    ]

    # in a number of datasets the iodine-containing molecules
    # were tainted due to an auxiliary basis set issue
    # This has since been resolved and entries have been recomputed
    # in new datasets, but we still need to filter the old ones
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
    if include_iodine:
        elements.append("I")

    # filter out other unsuitable entries
    dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        # HydrogenBondFilter(method="baker-hubbard"),
        ConnectivityFilter(tolerance=1.2),
        # UnperceivableStereoFilter(),
        # ElementFilter(allowed_elements=elements),
        ChargeCheckFilter(),
    )

    return dataset
