# Combine the sage 2.1 dataset with Pavan's

from openff.qcsubmit.results import (
    TorsionDriveResultCollection,
    OptimizationResultCollection,
)
from openff.toolkit import ForceField
from curate_dataset import select_parameters
import numpy as np
import json
import typing

Output: typing.TypeAlias = typing.Union[
    "TorsionDriveResultCollection", "OptimizationResultCollection"
]


def combine_datasets(d1: Output, d2: Output) -> Output:
    "Combine `d1` and `d2` and filter the resulting collection for duplicates"
    key = list(d1.entries.keys())[0]

    all_entries = d1.entries[key] + d2.entries[key]

    output_type = type(d1)

    assert output_type == type(d2)

    # filter duplicates
    unique = {record.record_id: record for record in all_entries}
    new_dataset = (output_type)(entries={key: list(unique.values())})

    return new_dataset


def combine_td(ff):
    base = "/home/brent/omsf/clone/sage-2.1.0/inputs-and-outputs/data-sets/"
    sage_td = TorsionDriveResultCollection.parse_file(
        base + "td-set-for-fitting-2.1.0.json"
    )

    records_to_remove = np.loadtxt("td_records_to_remove.dat", dtype=str)
    key = list(sage_td.entries.keys())[0]
    sage_td.entries[key] = [
        entry
        for entry in sage_td.entries[key]
        if entry.record_id not in records_to_remove
    ]

    pavan_td = TorsionDriveResultCollection.parse_file(
        "output/pavan-td-training-set.json"
    )

    combined_td = combine_datasets(sage_td, pavan_td)

    with open("output/combined-td.json", "w") as out:
        out.write(combined_td.json(indent=2))

    explicit_ring_torsions = np.loadtxt(
        "explicit_ring_torsions.dat", dtype=str
    )

    selected_parameters = select_parameters(
        combined_td,
        ["ProperTorsions"],
        force_field=ff,
        explicit_ring_torsions=explicit_ring_torsions,
    )
    with open("output/combined-td-smirks.json", "w") as file:
        json.dump(selected_parameters, file, indent=2)


def combine_opt(ff):
    base = "/home/brent/omsf/clone/sage-2.1.0/inputs-and-outputs/data-sets/"

    sage_td = OptimizationResultCollection.parse_file(
        base + "opt-set-for-fitting-2.1.0.json"
    )

    records_to_remove = np.loadtxt("opt_records_to_remove.dat", dtype=str)
    key = list(sage_td.entries.keys())[0]
    sage_td.entries[key] = [
        entry
        for entry in sage_td.entries[key]
        if entry.record_id not in records_to_remove
    ]

    pavan_td = OptimizationResultCollection.parse_file(
        "output/pavan-opt-training-set.json"
    )

    combined_td = combine_datasets(sage_td, pavan_td)

    with open("output/combined-opt.json", "w") as out:
        out.write(combined_td.json(indent=2))

    selected_parameters = select_parameters(
        combined_td,
        ["Bonds", "Angles"],
        force_field=ff,
    )
    with open("output/combined-opt-smirks.json", "w") as file:
        json.dump(selected_parameters, file, indent=2)


ff = ForceField(
    "../01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml",
    allow_cosmetic_attributes=True,
)

combine_td(ff)
combine_opt(ff)