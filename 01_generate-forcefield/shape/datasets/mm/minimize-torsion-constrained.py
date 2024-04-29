import qcportal
import pathlib

from collections import defaultdict
import typing

import click
from click_option_group import optgroup

import tqdm

import numpy as np

def get_best_all_atom_rmsd(molecule, qm_coordinates, mm_coordinates):
    from openff.toolkit.topology import Molecule
    from openff.units import unit
    from openeye import oechem

    molecule1 = Molecule(molecule)
    molecule1._conformers = [
        np.array(qm_coordinates).reshape((-1, 3)) * unit.angstrom
    ]
    molecule2 = Molecule(molecule)
    molecule2._conformers = [
        np.array(mm_coordinates).reshape((-1, 3)) * unit.angstrom
    ]
    return oechem.OERMSD(
        molecule1.to_openeye(),
        molecule2.to_openeye(),
        True,  # automorph
        False, # heavyOnly
        True,  # overlay
    )

def minimize_constrained(
    molecule,
    forcefield,
    dihedral = None,
):
    import openmm
    import openmm.unit
    from openff.units import unit
    from openff.units.openmm import from_openmm
    from openff.toolkit import Molecule, ForceField
    from openff.interchange.operations.minimize import _DEFAULT_ENERGY_MINIMIZATION_TOLERANCE

    forcefield = ForceField(forcefield, allow_cosmetic_attributes=True)
    interchange = forcefield.create_interchange(molecule.to_topology())

    simulation = interchange.to_openmm_simulation(
        openmm.LangevinMiddleIntegrator(
            293.15 * openmm.unit.kelvin,
            1.0 / openmm.unit.picosecond,
            2.0 * openmm.unit.femtosecond,
        ),
        combine_nonbonded_forces=True,
    )

    simulation.context.computeVirtualSites()

    if dihedral is not None:
        for i in dihedral:
            simulation.system.setParticleMass(i, 0.0)

    simulation.minimizeEnergy(
        tolerance=_DEFAULT_ENERGY_MINIMIZATION_TOLERANCE.to_openmm(),
        maxIterations=10_000,
    )

    coordinates = from_openmm(
        simulation.context.getState(getPositions=True).getPositions(asNumpy=True)[
            : interchange.positions.shape[0],
            :,
        ],
    )
    return coordinates



def benchmark_single(
    row,
    force_field: str,
):
    from openff.toolkit import Molecule, ForceField
    from openff.units import unit
    from openff.interchange.drivers.openmm import get_openmm_energies
    from ibstore.analysis import get_rmsd, get_tfd, get_internal_coordinate_rmsds

    molecule = Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True,
    )
    molecule._conformers = [
        np.array(row["conformer"]).reshape((-1, 3)) * unit.angstrom
    ]

    ff = ForceField(force_field, allow_cosmetic_attributes=True)
    parameter_labels = {}
    for parameter_name, handler in ff._parameter_handlers.items():
        for parameter in handler.parameters:
            key = f"{parameter_name}_{parameter.smirks}"
            parameter_labels[key] = False
    
    labels = ff.label_molecules(molecule.to_topology())[0]
    for parameter_name, assigned in labels.items():
        for value in assigned.values():
            key = f"{parameter_name}_{value.smirks}"
            parameter_labels[key] = True
    for pattern in ["[r4:1]", "[#7]-[#16:1](=[#8])(=[#8])"]:
        parameter_labels[pattern] = bool(molecule.chemical_environment_matches(pattern))

    minimized_positions = minimize_constrained(
        molecule,
        force_field,
        dihedral=row["dihedral"],
    )

    molecule._conformers = [minimized_positions]
    ic2 = ff.create_interchange(molecule.to_topology())
    energies = get_openmm_energies(ic2)


    n_heavy_atoms = sum(1 for atom in molecule.atoms if atom.atomic_number > 1)

    entry = {
        "forcefield": force_field,
        "qcarchive_id": row["qcarchive_id"],
        "mapped_smiles": row["mapped_smiles"],
        "n_atoms": len(molecule.atoms),
        "n_heavy_atoms": n_heavy_atoms,
        "mm_coordinates": minimized_positions.m_as(unit.angstrom).flatten().tolist(),
        "qm_coordinates": row["conformer"],
        "qm_energy": row["energy"],
        "mm_energy": energies.total_energy.m_as(unit.kilocalories_per_mole),
    }
    for key, value in energies.energies.items():
        entry[f"energy_{key}"] = value.m_as(unit.kilocalories_per_mole)

    molecule =  Molecule.from_mapped_smiles(
        row["mapped_smiles"],
        allow_undefined_stereo=True,
    )
    qm_coordinates = np.array(row["conformer"]).reshape((-1, 3))
    mm_coordinates = minimized_positions.m_as(unit.angstrom).reshape((-1, 3))
    entry["RMSD"] = get_rmsd(molecule, qm_coordinates, mm_coordinates)
    entry["RMSD_AA"] = get_best_all_atom_rmsd(molecule, qm_coordinates, mm_coordinates)
    try:
        entry["TFD"] = get_tfd(molecule, qm_coordinates, mm_coordinates)
    except IndexError:
        entry["TFD"] = np.nan
    internal_coordinate_rmsds = get_internal_coordinate_rmsds(
        molecule, qm_coordinates, mm_coordinates
    )
    entry.update(internal_coordinate_rmsds)
    entry.update(parameter_labels)



    return entry


def batch_compute(
    rows: list[dict],
    force_field: str = None,
):

    results = []
    errors = []

    for row in tqdm.tqdm(rows):
        try:
            inchi_results = benchmark_single(row, force_field)
        except Exception as e:
            errors.append(f"{row['torsiondrive_id']}: {row['qcarchive_id']} --- {str(e)}")
            # raise e
        else:
            results.append(inchi_results)
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
    input_dataset_path,
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

    dataset = ds.dataset(input_dataset_path, format="parquet")
    columns = [
        "mapped_smiles",
        "conformer",
        "qcarchive_id",
        "torsiondrive_id", 
        "energy",
        "dihedral",
    ]
    entries = dataset.to_table(columns=columns).to_pylist()

    output_directory = pathlib.Path(output_dataset)
    output_directory.mkdir(exist_ok=True, parents=True)
    start_index = 0
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


    all_data_entries = []
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
            table_path = output_directory / f"batch-{i:04d}.parquet"

            pq.write_table(batch_table, table_path)
            print(f"Wrote {len(batch)} to {table_path}")
            
            all_errors.extend(errors)
            all_data_entries.extend(batch)

    with open("errors.txt", "w") as file:
        file.write("\n".join(all_errors))


if __name__ == "__main__":
    compute_benchmarks()
