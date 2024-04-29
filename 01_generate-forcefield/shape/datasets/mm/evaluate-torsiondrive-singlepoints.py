import qcportal
from collections import defaultdict
import typing

import pathlib
import click
from click_option_group import optgroup

import tqdm

import numpy as np
from openff.toolkit import Molecule, ForceField

def benchmark_single(
    row,
    force_field: ForceField,
):
    from openff.units import unit
    from openff.toolkit import Molecule, ForceField
    from openff.interchange.drivers import get_openmm_energies

    entry = dict(row)
    entry["qm_energy"] = entry.pop("energy")
    entry["qm_coordinates"] = entry["conformer"]

    molecule = Molecule.from_mapped_smiles(
        entry["mapped_smiles"],
        allow_undefined_stereo=True,
    )
    molecule._conformers = [
        np.array(entry["conformer"]).reshape((-1, 3)) * unit.angstrom
    ]
    ic = force_field.create_interchange(molecule.to_topology())
    energies = get_openmm_energies(ic, detailed=True, combine_nonbonded_forces=False)
    for k, v in energies.energies.items():
        entry[k] = v.m_as(unit.kilocalories_per_mole)
    mm_energy = energies.total_energy.m_as(unit.kilocalories_per_mole)
    entry["mm_energy"] = mm_energy
    return entry


def batch_compute(
    rows: list[dict],
    force_field: str = None,
):

    results = []
    errors = []

    ff = ForceField(force_field, allow_cosmetic_attributes=True)
    for row in tqdm.tqdm(rows):
        try:
            result = benchmark_single(row, ff)
        except Exception as e:
            errors.append(str(e))
        else:
            result["forcefield"] = force_field
            results.append(result)
    return results, errors

@click.command()
@click.option(
    "--input",
    "input_dataset_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help="The path to the input Parquet path of the torsiondrive dataset. ",
    default="qm/output"
)
@click.option(
    "--force-field",
    type=str,
    help="The name or path of the force field to use to optimize in the database.",
)
@click.option(
    "--output",
    "output_dataset",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    help="The path to the output Parquet path.",
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("Cluster configuration", help="Options to configure cluster workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The SLURM queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that SLURM workers should run using.",
    type=str,
)
def compute_benchmarks(
    force_field: str,
    output_dataset: str,
    input_dataset_path: str = None,
    worker_type: typing.Literal["slurm", "local"] = "local",
    queue: str = "free",
    conda_environment: str = "ib-dev",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    import pyarrow.dataset as ds

    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed

    output_directory = pathlib.Path(output_dataset)
    output_directory.mkdir(exist_ok=True, parents=True)
    start_index = 0

    dataset = ds.dataset(input_dataset_path, format="parquet")
    entries = dataset.to_table().to_pylist()

    try:
        print(f"Querying {output_directory}")
        existing = ds.dataset(output_directory)
        if existing.count_rows():
            expression = pc.field("forcefield") == force_field
            subset = existing.filter(expression)
            qca_ids = subset.to_table(columns=["qcarchive_id"]).to_pydict()["qcarchive_id"]
            entries = [
                entry for entry in entries if entry["qcarchive_id"] not in qca_ids
            ]

            start_index = len(existing.files)
    except BaseException as e:
        print(e)
    
    print(f"Starting from start_index={start_index}")

    all_errors = []

    with batch_distributed(
        entries,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(
            batch_compute,
            force_field=force_field,
        ))
        futures = list(batcher(
            batch_compute,
            force_field=force_field,
        ))
        for i, future in tqdm.tqdm(
            enumerate(
                distributed.as_completed(futures, raise_errors=False),
                start_index,
            ),
            total=len(futures),
            desc="Updating entries",
        ):
            batch, errors = future.result()

            batch_table = pa.Table.from_pylist(batch)
            print(batch_table.schema)
            table_path = output_directory / f"batch-{i:04d}.parquet"

            pq.write_table(batch_table, table_path)
            print(f"Wrote {len(batch)} to {table_path}")
            
            all_errors.extend(errors)

    with open(output_directory / "errors.txt", "w") as file:
        file.write("\n".join(all_errors))
 


if __name__ == "__main__":
    compute_benchmarks()
