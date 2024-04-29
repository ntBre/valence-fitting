# Datasets

This was generated outside this repo on HPC3.

The scripts here download QCArchive datasets and convert them to Parquet format.

## Files

- run-*.sh files show intended usage of Python scripts.
- conda_env.yaml shows the full Python environment used to run the Python scripts
- output/torsiondrive/* contains output data.

## Schema

* `dataset` (string): Dataset name of record
* `qcarchive_id` (int64): QCArchive ID of particular geometry
* `torsiondrive_id` (int64): ID of particular TorsionDrive (which may contain many geometries)
* `mapped_smiles` (string): The CMILES attached to the QCArchive dataset molecule. This is *not* cleaned by roundtripping through the Toolkit, but is read from the data
* `n_atoms` (int64): Number of atoms in the molecule
* `energy` (double): Energy of the geometry, in kcal/mol
* `conformer` (list<double>): a flattened list of the geometry in angstrom. This should be reshaped to (n_atoms, 3) for use
* `atomic_numbers` (list<int64>): a list of unique atomic numbers in the molecule
* `grid_id` (int64): the value of the torsion angle in degrees
* `dihedral` (list<int64>): the atom indices (indexed from 0) of the torsion in the mapped smiles
* `atomic_number_1` (<int64>): the atomic number of the *second* atom in `dihedral`, i.e. one of the central bond atoms
* `atomic_number_2` (<int64>): the atomic number of the *third* atom in `dihedral`, i.e. one of the central bond atoms
