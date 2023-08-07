from collections import defaultdict, Counter
import functools
import json
import logging
import random
import typing
import os.path
import sys

import numpy as np
import click

from qcportal.models import TorsionDriveRecord, OptimizationRecord
from openff.qcsubmit.results import (
    TorsionDriveResultCollection,
    OptimizationResultCollection,
)
from openff.toolkit import ForceField, Molecule


def check_torsion_is_in_ring(
    molecule: "Molecule",
    indices: typing.Tuple[int, int, int, int],
) -> bool:
    """
    Check if a torsion is in a ring.

    If a torsion I-J-K-L is given, it checks
    whether all bonds I-J, J-K, and K-L are in a ring.

    """
    i, j, k, m = indices
    return (
        molecule.get_bond_between(i, j).is_in_ring()
        and molecule.get_bond_between(j, k).is_in_ring()
        and molecule.get_bond_between(k, m).is_in_ring()
    )


def label_and_tag_ids(
    record_and_molecule: typing.Tuple[
        typing.Union["TorsionDriveRecord", "OptimizationRecord"], "Molecule"
    ],
    force_field: "ForceField",
    parameter_types: typing.List[str],
    explicit_ring_torsions: typing.Optional[str] = None,
) -> typing.Set[typing.Tuple[str, str, int]]:
    from qcportal.models import TorsionDriveRecord

    if explicit_ring_torsions is not None:
        ring_torsions = np.loadtxt(explicit_ring_torsions, dtype=str)
    else:
        ring_torsions = []

    record, molecule = record_and_molecule
    mol_labels = force_field.label_molecules(molecule.to_topology())[0]
    parameter_ids = set()

    for parameter_type in parameter_types:
        parameter_labels = mol_labels[parameter_type]

        for indices, parameter in parameter_labels.items():
            # remove mismatching torsiondrives
            if isinstance(record, TorsionDriveRecord):
                # check central bond, i.e. middle 2 atoms
                record_atoms = record.keywords.dihedrals[0]
                if set(indices[1:3]) != set(record_atoms[1:3]):
                    continue

                # some general parameters overlap with in-ring torsions and
                # there are many torsion scans from Gen1 sets that have
                # in-ring torsions and we want to exclude them in training
                # as they result in higher k values unless the parameters
                # have smirks explicitly for an in-ring torsion. It is to be
                # noted that training on in-ring torsions is needed to
                # properly model puckering in rings with hetero atoms
                if parameter.id not in ring_torsions:
                    if check_torsion_is_in_ring(molecule, indices):
                        continue

            n_heavy_atoms = sum(
                1 for atom in molecule.atoms if atom.atomic_number != 1
            )
            parameter_ids.add((parameter.id, record.id, n_heavy_atoms))
    return parameter_ids


def get_parameter_distribution(
    dataset: typing.Union[
        "TorsionDriveResultCollection", "OptimizationResultCollection"
    ],
    parameter_types: typing.List[str],
    force_field: "ForceField",
    explicit_ring_torsions: typing.Optional[str] = None,
    n_processes: int = 4,
) -> typing.Tuple[
    Counter, typing.Dict[str, typing.List[typing.Tuple[int, str]]]
]:
    coverage = Counter()
    parameter_records = defaultdict(list)

    func = functools.partial(
        label_and_tag_ids,
        force_field=force_field,
        parameter_types=parameter_types,
        explicit_ring_torsions=explicit_ring_torsions,
    )
    for record in dataset.to_records():
        parameter_ids = func(record)
        for parameter_id, record_id, n_heavy_atoms in parameter_ids:
            coverage[parameter_id] += 1
            parameter_records[parameter_id].append((n_heavy_atoms, record_id))

    return coverage, dict(parameter_records)


def cap_torsions_per_parameter(
    force_field: "ForceField",
    dataset: "TorsionDriveResultCollection",
    cap_size: int = 5,
    explicit_ring_torsions: typing.Optional[str] = None,
    method: typing.Literal[
        "pick_random", "pick_heavy", "pick_light"
    ] = "pick_random",
    verbose: bool = True,
    n_processes: int = 4,
):
    coverage, parameter_records = get_parameter_distribution(
        dataset=dataset,
        parameter_types=["ProperTorsions"],
        force_field=force_field,
        explicit_ring_torsions=explicit_ring_torsions,
        n_processes=n_processes,
    )
    records_to_keep = {}
    for parameter_id in coverage:
        if coverage[parameter_id] <= cap_size:
            n_atom_records = parameter_records[parameter_id]
        else:
            if method == "pick_heavy":
                n_atom_records = sorted(
                    parameter_records[parameter_id],
                    key=lambda x: x[0],
                    reverse=True,
                )[:cap_size]
            elif method == "pick_light":
                n_atom_records = sorted(
                    parameter_records[parameter_id],
                    key=lambda x: x[0],
                    reverse=False,
                )[:cap_size]
            elif method == "pick_random":
                n_atom_records = random.sample(
                    parameter_records[parameter_id], cap_size
                )

        _, records = zip(*n_atom_records)
        records_to_keep[parameter_id] = records

    if verbose:
        print("Final coverage")
        for parameter_id, records in records_to_keep.items():
            print(
                f"{parameter_id:>6s}: {len(records):>4d} "
                f"/ {coverage[parameter_id]:>4d} records"
            )

    ids_to_keep = [
        record_id
        for record_ids in records_to_keep.values()
        for record_id in record_ids
    ]
    print(f"Total records: {dataset.n_results}")
    print(f"Total records to keep: {len(ids_to_keep)}")

    key = list(dataset.entries.keys())[0]
    dataset.entries[key] = [
        record
        for record in dataset.entries[key]
        if record.record_id in ids_to_keep
    ]
    return dataset


def download_td_data(
    td_datasets: typing.List[str],
    ds_cache,
    invalidate_cache=False,
) -> "TorsionDriveResultCollection":
    """Download TorsionDrive datasets."""

    from qcportal import FractalClient
    from openff.qcsubmit.results import TorsionDriveResultCollection

    if os.path.isfile(ds_cache) and not invalidate_cache:
        print(f"loading td from {ds_cache}", file=sys.stderr)
        return TorsionDriveResultCollection.parse_file(ds_cache)

    client = FractalClient()
    dataset = TorsionDriveResultCollection.from_server(
        client=client,
        datasets=td_datasets,
        spec_name="default",
    )
    with open(ds_cache, "w") as out:
        out.write(dataset.json(indent=2))

    return dataset


def select_parameters(
    dataset: typing.Union[
        "TorsionDriveResultCollection", "OptimizationResultCollection"
    ],
    parameter_types: typing.List[str],
    force_field: "ForceField",
    explicit_ring_torsions: typing.Optional[str] = None,
    n_processes: int = 1,
    min_coverage: int = 5,
):
    coverage, _ = get_parameter_distribution(
        dataset=dataset,
        parameter_types=parameter_types,
        force_field=force_field,
        explicit_ring_torsions=explicit_ring_torsions,
        n_processes=n_processes,
    )

    selected_parameters = defaultdict(list)
    for parameter_type in parameter_types:
        handler = force_field.get_parameter_handler(parameter_type)

        for parameter_id, count in coverage.items():
            if count < min_coverage:
                continue
            parameters = handler.get_parameter({"id": parameter_id})
            if not len(parameters):
                continue
            selected_parameters[parameter_type].append(parameters[0].smirks)
    return selected_parameters


@click.group()
def cli():
    pass


@cli.command("download-td")
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--output-parameter-smirks",
    "output_parameter_smirks_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--initial-forcefield",
    required=True,
    type=str,
    help=(
        "The name of the initial force field to use. "
        "Alternatively, the path to a force field"
    ),
)
@click.option(
    "--explicit-ring-torsions",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help=(
        "The path to a file containing a list of parameter IDs that are "
        "ring torsions. This should be a text file with one ID per line."
    ),
)
@click.option(
    "--cap-size",
    type=int,
    default=5,
    show_default=True,
    help=(
        "The maximum number of torsions to include per parameter "
        "in the auxiliary datasets."
        "If there are more torsions than this, a subset will be selected."
    ),
)
@click.option(
    "--cap-method",
    type=click.Choice(["pick_random", "pick_heavy", "pick_light"]),
    default="pick_random",
    show_default=True,
    help=(
        "The method to use to select the torsions to include per parameter "
        "in the auxiliary datasets."
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to print out additional information.",
)
@click.option(
    "--n-processes",
    type=int,
    default=4,
    show_default=True,
    help="The number of processes to use when processing the data.",
)
@click.option(
    "--min-record-coverage",
    type=int,
    default=5,
    show_default=True,
    help=(
        "The minimum number of records a parameter must have to be included "
        "in the force field optimization."
    ),
)
def get_td_data(
    output_path: str,
    output_parameter_smirks_path: str,
    initial_forcefield: str,
    explicit_ring_torsions: typing.Optional[str] = None,
    cap_size: int = 5,
    cap_method: typing.Literal[
        "pick_random", "pick_heavy", "pick_light"
    ] = "pick_random",
    verbose: bool = True,
    n_processes: int = 4,
    min_record_coverage: int = 5,
):
    """
    Download TorsionDrive data in stages.
    \f

    1. Download the core datasets and filter out unsuitable entries.
    2. Download the auxiliary datasets and filter out unsuitable entries.
    3. Cap the number of auxiliary torsions per parameter to a maximum
       of ``cap_size``. This can be done by picking random torsions, or
       selecting those with the least (``pick_light``) or most
       (``pick_heavy``) heavy atoms.
    4. Add additional torsiondrive records from a file.
    5. Filter out duplicate torsiondrive records.
    6. Filter out molecules that fail AM1-BCC ELF10 charging. This step
       is the slowest so it's done last.

    The unsuitability filters are:
        - removing incomplete entries
        - removing entries with hydrogen bonds
        - removing entries with unperceivable stereochemistry
        - removing entries with connectivity rearrangements
        - removing entries with iodine

    Parameters
    ----------
    core_td_datasets
        The core torsiondrive datasets to download.
        These are filtered for unsuitability, but are not capped.
    aux_td_datasets
        The auxiliary torsiondrive datasets to download.
        These are filtered for unsuitability, and are capped
        to a certain number of torsions per parameter.
    initial_forcefield
        The initial forcefield to use for filtering torsiondrive entries.
    td_records_to_remove
        A file containing a list of torsiondrive record IDs to remove.
    additional_td_records
        A file containing a list of additional torsiondrive records to add.
        This should be a JSON file of a ``TorsionDriveResultCollection``.
    cap_size
        The maximum number of torsions to keep per parameter.
    cap_method
        The method to use to cap the number of torsions per parameter.
        One of ``pick_random``, ``pick_heavy``, or ``pick_light``.
    verbose
        Whether to print out information about the number of records
        at each stage.
    n_processes
        The number of processes to use for multiprocessing.
    """
    from openff.toolkit import ForceField
    from openff.qcsubmit.results import TorsionDriveResultCollection

    # suppress stereochemistry warnings
    logging.getLogger("openff").setLevel(logging.ERROR)

    ff = ForceField(initial_forcefield, allow_cosmetic_attributes=True)

    core_dataset = TorsionDriveResultCollection.parse_file(
        "datasets/filtered-core-td.json"
    )

    key = list(core_dataset.entries.keys())[0]

    all_entries = core_dataset.entries[key]

    # filter in case we have doubled up records
    unique_entries = {record.record_id: record for record in all_entries}
    new_dataset = TorsionDriveResultCollection(
        entries={key: list(unique_entries.values())}
    )

    n = new_dataset.n_results
    print(f"final number of core td entries: {n}")

    with open(output_path, "w") as file:
        file.write(new_dataset.json(indent=2))
    if verbose:
        print(f"Saved to {output_path}")

    selected_parameters = select_parameters(
        new_dataset,
        ["ProperTorsions"],
        force_field=ff,
        explicit_ring_torsions=explicit_ring_torsions,
        n_processes=n_processes,
        min_coverage=min_record_coverage,
    )
    with open(output_parameter_smirks_path, "w") as file:
        json.dump(selected_parameters, file, indent=2)


@cli.command("download-opt")
@click.option(
    "--input-dataset",
    "input_dataset",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--output-parameter-smirks",
    "output_parameter_smirks_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--initial-forcefield",
    required=True,
    type=str,
    help=(
        "The name of the initial force field to use. "
        "Alternatively, the path to a force field"
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to print out additional information.",
)
@click.option(
    "--n-processes",
    type=int,
    default=4,
    show_default=True,
    help="The number of processes to use when processing the data.",
)
@click.option(
    "--min-record-coverage",
    type=int,
    default=5,
    show_default=True,
    help=(
        "The minimum number of records a parameter must have to be included "
        "in the force field optimization."
    ),
)
def get_opt_data(
    output_path: str,
    input_dataset: str,
    output_parameter_smirks_path: str,
    initial_forcefield: str,
    verbose: bool = True,
    n_processes: int = 4,
    min_record_coverage: int = 5,
):
    """Download and filter optimization datasets.

    \f
    Parameters
    ----------
    core_opt_datasets
        The core optimization datasets to download.
    opt_records_to_remove
        A file containing a list of optimization record IDs to remove.
    max_opt_conformers
        The maximum number of conformers to keep per molecule.
        Conformers are filled using a greedy RMSD filter
    """
    from openff.qcsubmit.results import OptimizationResultCollection
    from openff.toolkit import ForceField

    # suppress stereochemistry warnings
    logging.getLogger("openff").setLevel(logging.ERROR)

    ff = ForceField(initial_forcefield, allow_cosmetic_attributes=True)

    core_dataset = OptimizationResultCollection.parse_file(input_dataset)

    key = list(core_dataset.entries.keys())[0]

    # filter in case we have doubled up records
    unique_entries = {
        record.record_id: record for record in core_dataset.entries[key]
    }
    new_dataset = OptimizationResultCollection(
        entries={key: list(unique_entries.values())}
    )

    n = new_dataset.n_results
    print(f"final number of core opt entries: {n}")

    with open(output_path, "w") as file:
        file.write(new_dataset.json(indent=2))
    if verbose:
        print(f"Saved to {output_path}")

    selected_parameters = select_parameters(
        new_dataset,
        ["Bonds", "Angles"],
        force_field=ff,
        n_processes=n_processes,
        min_coverage=min_record_coverage,
    )
    with open(output_parameter_smirks_path, "w") as file:
        json.dump(selected_parameters, file, indent=2)


if __name__ == "__main__":
    cli()
