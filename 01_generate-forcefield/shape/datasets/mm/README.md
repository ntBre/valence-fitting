# Datasets

This was generated outside this repo on HPC3.

The scripts here minimize TorsionDrive conformers and evaluate the energies.

## Files

- run-*.sh files show intended usage of Python scripts.
- conda_env.yaml shows the full Python environment used to run the Python scripts
- output data is in `singlepoint-torsiondrive-datasets` and `minimized-torsiondrive-datasets`
- Note: some *.parquet files may be empty due to the way the scripts write intermediate files.

## Schema

### singlepoint-torsiondrive-datasets
- `dataset` (string): Dataset name of record
- `qcarchive_id` (int64): QCArchive ID of particular geometry
- `torsiondrive_id` (int64): ID of particular TorsionDrive (which may contain many geometries)
- `mapped_smiles` (string): The CMILES attached to the QCArchive dataset molecule. This is *not* cleaned by roundtripping through the Toolkit, but is read from the data
- `n_atoms` (int64): Number of atoms in the molecule
- `atomic_numbers` (list<int64>): a list of unique atomic numbers in the molecule
- `grid_id` (int64): the value of the torsion angle in degrees
- `dihedral` (list<int64>): the atom indices (indexed from 0) of the torsion in the mapped smiles
- `atomic_number_1` (<int64>): the atomic number of the *second* atom in `dihedral`, i.e. one of the central bond atoms
- `atomic_number_2` (<int64>): the atomic number of the *third* atom in `dihedral`, i.e. one of the central bond atoms
- `N` (bool): whether this molecule has this element
- `O` (bool): whether this molecule has this element
- `F` (bool): whether this molecule has this element
- `P` (bool): whether this molecule has this element
- `S` (bool): whether this molecule has this element
- `Cl` (bool): whether this molecule has this element
- `Br` (bool): whether this molecule has this element
- `I` (bool): whether this molecule has this element
- `qm_energy` (double): QM energy of the geometry, in kcal/mol
- `qm_coordinates` (list<double>): a flattened list of the geometry in angstrom. This should be reshaped to (n_atoms, 3) for use
- `vdW` (double): Component of energy in `mm_energy`, in kcal/mol
- `Electrostatics` (double): Component of energy in `mm_energy`, in kcal/mol
- `vdW 1-4` (double): Component of energy in `mm_energy`, in kcal/mol
- `Electrostatics 1-4` (double): Component of energy in `mm_energy`, in kcal/mol
- `Torsion` (double): Component of energy in `mm_energy`, in kcal/mol
- `Angle` (double): Component of energy in `mm_energy`, in kcal/mol
- `Bond` (double): Component of energy in `mm_energy`, in kcal/mol
- `mm_energy` (double): MM energy of the geometry, in kcal/mol
- `forcefield` (string): Force field used to evaluate


### minimized-torsiondrive-datasets

- `forcefield` (string): Force field used to minimize and evaluate
- `qcarchive_id` (int64): QCArchive ID of particular QM geometry
- `mapped_smiles` (string): The CMILES attached to the QCArchive data
- `n_atoms` (int64):  Number of atoms in the molecule
- `n_heavy_atoms` (int64): Number of heavy atoms in the molecule
- `mm_coordinates` (list<double>): a flattened list of the geometry in angstrom. This should be reshaped to (n_atoms, 3) for use
- `qm_coordinates` (list<double>): a flattened list of the geometry in angstrom. This should be reshaped to (n_atoms, 3) for use
- `qm_energy` (double): QM energy of QM coordinates in kcal/mol
- `mm_energy` (double): MM energy of MM coordinates in kcal/mol
- `energy_Bond` (double): Bond component of mm_energy in kcal/mol
- `energy_Angle` (double): Angle component of mm_energy in kcal/mol
- `energy_Torsion` (double): Torsion component of mm_energy in kcal/mol
- `energy_Nonbonded` (double): Nonbonded component of mm_energy in kcal/mol
- `RMSD` (double): heavy-atom RMSD (from YAMMBS)
- `RMSD_AA` (double): all-atom RMSD
- `TFD` (double): TFD (from YAMMBS), if any
- `Bond` (double): Internal coordinate RMSD (from YAMMBS)
- `Angle` (double): Internal cooordinate RMSD (from YAMMBS)
- `Dihedral` (double): Internal cooordinate RMSD (from YAMMBS)
- `Bonds_*` (bool): presence of this particular Bond parameter in molecule
- `Angles_*` (bool): presence of this particular Angle parameter in molecule
- ... etc for all parameters in the `forcefield`
