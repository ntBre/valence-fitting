import logging
import typing
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusEnum,
    RecordStatusFilter,
    SinglepointRecordFilter,
    T,
)
from openff.toolkit.utils.exceptions import (
    ChargeCalculationError,
    ConformerGenerationError,
)
from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def imap_fn(record_and_molecule):
    """record here is actually only a record.id because the records maintain a
    reference to a PortalClient, which cannot be pickled for multiprocessing.
    """
    record, molecule = record_and_molecule
    try:
        OpenEyeToolkitWrapper().assign_partial_charges(
            molecule, partial_charge_method="am1bccelf10"
        )
    except (ChargeCalculationError, ConformerGenerationError):
        ok = False
    else:
        ok = True

    return record, ok


class ChargeCheckFilter(SinglepointRecordFilter):
    nprocs: int = 1
    chunksize: int = 1

    # this is not needed now that I overrode _apply, and I need to pass along
    # the record_id anyway, but pydantic requires it to be here
    def _filter_function(self, result, record, molecule) -> bool:
        raise NotImplementedError()

    def _apply(self, result_collection: T) -> T:
        """Copy of SinglepointRecordFilter._apply with added logging, progress
        reporting, and eventually parallelism."""

        all_records_and_molecules = defaultdict(list)

        logger.info("starting to_records")

        for record, molecule in result_collection.to_records():
            all_records_and_molecules[record._client.address].append(
                (record.id, molecule)
            )

        logger.info("finished to records")

        filtered_results = {}

        for address, entries in result_collection.entries.items():
            records_and_molecules = all_records_and_molecules[address]

            filtered_ids = []
            with Pool(processes=self.nprocs) as p:
                for record_id, ok in tqdm(
                    p.imap(
                        imap_fn,
                        records_and_molecules,
                        chunksize=self.chunksize,
                    ),
                    total=len(records_and_molecules),
                    desc="Filtering charge errors",
                ):
                    if ok:
                        filtered_ids.append(record_id)

            filtered_results[address] = [
                entry for entry in entries if entry.record_id in filtered_ids
            ]

        result_collection.entries = filtered_results

        return result_collection


class NoisyFilter(SinglepointRecordFilter):
    """A filter that always returns true but can be used for signaling progress
    through a sequence of real filters."""

    name: str

    def _apply(self, result_collection):
        n = result_collection.n_results
        print(f"starting filter: {self.name} on {n} records")
        return super()._apply(result_collection)

    def _filter_function(self, result, record, molecule) -> bool:
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
