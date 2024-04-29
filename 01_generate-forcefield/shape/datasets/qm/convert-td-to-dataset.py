import qcportal # necessary for zstd mismatch error
import json

import click
import tqdm
import qcelemental

HARTREE2KCALMOL = qcelemental.constants.hartree2kcalmol
BOHR2ANGSTROM = qcelemental.constants.bohr2angstroms

@click.command()
@click.option(
    "--input",
    "input_file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="The path to the input JSON file of OptimizationRecords.",
)
@click.option(
    "--name",
    "dataset_name",
    type=str,
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    help="The path to the output Parquet path of the dataset.",
)
def convert(
    input_file: str,
    dataset_name: str,
    output_path: str,
):
    from openff.toolkit import Molecule
    from openff.qcsubmit.results import TorsionDriveResultCollection

    import pyarrow as pa
    import pyarrow.dataset as ds

    dataset = TorsionDriveResultCollection.parse_file(input_file)
    td_record_to_cmiles = {
        entry.record_id: entry.cmiles
        for entry in dataset.entries['https://api.qcarchive.molssi.org:443/']
    }

    ELEMENTS = {
        "N": 7,
        "O": 8,
        "F": 9,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Br": 35,
        "I": 53
    }

    records_and_molecules = dataset.to_records()
    all_entries = []
    for record, openff_molecule in tqdm.tqdm(records_and_molecules):
        cmiles = td_record_to_cmiles[record.id]
        all_atomic_numbers = [
            atom.atomic_number
            for atom in openff_molecule.atoms
        ]
        atomic_numbers = sorted(set(all_atomic_numbers))
        els = {
            k: v in atomic_numbers
            for k, v in ELEMENTS.items()
        }
        dihedrals = record.specification.keywords.dihedrals[0]

        central_bond = sorted(
            [
                all_atomic_numbers[dihedrals[1]],
                all_atomic_numbers[dihedrals[2]]
            ]
        )
        
        for grid_id, optimization in record.minimum_optimizations.items():
            conformer = optimization.final_molecule.geometry * BOHR2ANGSTROM
            entry = {
                "dataset": dataset_name,
                "qcarchive_id": optimization.id,
                "torsiondrive_id": record.id,
                "mapped_smiles": cmiles,
                "n_atoms": len(openff_molecule.atoms),
                "energy": optimization.energies[-1] * HARTREE2KCALMOL,
                "conformer": conformer.flatten().tolist(),
                "atomic_numbers": atomic_numbers,
                "grid_id": grid_id[0],
                "dihedral": list(dihedrals),
                "atomic_number_1": central_bond[0],
                "atomic_number_2": central_bond[1]
            }
            entry.update(els)
            all_entries.append(entry)
    
    table = pa.Table.from_pylist(all_entries)
    print(table.schema)
    ds.write_dataset(table, output_path, format="parquet")
    print(f"Wrote {len(all_entries)} to {output_path}")



if __name__ == "__main__":
    convert()