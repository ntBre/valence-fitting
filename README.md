# Valence fitting

This repository is a small modification of scripts at https://github.com/openforcefield/sage-2.1.0.

Create the environment necessary with the below. Replace ``micromamba`` with ``mamba`` or ``conda`` as needed.

```bash
micromamba env create -f environment.yaml
```

If you're running on Apple Silicon, some packages are still not available for ARM yet. Instead, do the following:


```bash
CONDA_SUBDIR=osx-64 micromamba env create -f environment.yaml
micromamba activate fb-195-tk-014-py310
conda config --env --set subdir osx-64  # make sure updates keep using the subdir
```

The general steps laid out here are:

1. Generate the initial force field for fitting.
   This step defines the SMIRKS of each parameter that will be fit.
   When adding, removing, or modifying parameters in force fields, here
   is where it should be changed.
2. Download and curate QM datasets.
   This step downloads QM datasets; manually adds additional data records;
   and manually removes bad records. For torsions, you can specify whether
   to cap the number of torsion data records per torsion parameter for fitting.
   For optimizations, you can cap the number of conformers used as well.
   The output of these scripts are the QM datasets, as well as which parameters
   will be trained. A parameter will not be trained if there isn't
   enough data coverage (also specified by the user).
3. This step generates initial values for each force field parameter
   using the Modified Seminario Method.
4. This step sets up the final ForceBalance input files for training.

After setting up ForceBalance, the ``scripts`` directory contains scripts
suitable for running the fit on HPC3.


Each of the scripts in the numbered directories can be used via commandline,
as demonstrated in the accompanying bash scripts:

```bash
python curate-dataset.py --help
```
```
Usage: curate-dataset.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  download-opt  Download and filter optimization datasets.
  download-td   Download TorsionDrive data in stages.
```
```bash
python curate-dataset.py download-opt --help
```
```
Usage: curate-dataset.py download-opt [OPTIONS]

  Download and filter optimization datasets.

Options:
  --output FILE                   The path to write the dataset to. Should be
                                  a JSON  [required]
  --output-parameter-smirks FILE  The path to write the dataset to. Should be
                                  a JSON  [required]
  --initial-forcefield TEXT       The name of the initial force field to use.
                                  Alternatively, the path to a force field
                                  [required]
  --core-opt-dataset TEXT         The name of an optimization dataset to
                                  download. These will have iodine molecules
                                  filtered out.  [required]
  --iodine-opt-dataset TEXT       The name of an optimization dataset to
                                  download. These will have iodine molecules
                                  included.  [required]
  --opt-records-to-remove FILE    The path to a file containing a list of
                                  record IDs to remove. This should be a text
                                  file with one record ID per line.
  --verbose                       Whether to print out additional information.
  --max-opt-conformers INTEGER    The maximum number of conformers to keep per
                                  molecule.  [default: 12]
  --n-processes INTEGER           The number of processes to use when
                                  processing the data.  [default: 4]
  --min-record-coverage INTEGER   The minimum number of records a parameter
                                  must have to be included in the force field
                                  optimization.  [default: 5]
  --help                          Show this message and exit.
```
```bash
python curate-dataset.py download-td --help
```
```
Usage: curate-dataset.py download-td [OPTIONS]

  Download TorsionDrive data in stages.

Options:
  --output FILE                   The path to write the dataset to. Should be
                                  a JSON  [required]
  --output-parameter-smirks FILE  The path to write the dataset to. Should be
                                  a JSON  [required]
  --core-td-dataset TEXT          The name of a torsiondrive dataset to
                                  download.  [required]
  --aux-td-dataset TEXT           The name of a torsiondrive dataset to
                                  download.  [required]
  --initial-forcefield TEXT       The name of the initial force field to use.
                                  Alternatively, the path to a force field
                                  [required]
  --explicit-ring-torsions FILE   The path to a file containing a list of
                                  parameter IDs that are ring torsions. This
                                  should be a text file with one ID per line.
  --td-records-to-remove FILE     The path to a file containing a list of
                                  record IDs to remove. This should be a text
                                  file with one record ID per line.
  --additional-td-records FILE    The path to a file containing a
                                  TorsionDriveResultCollection containing
                                  additional torsiondrive records to include.
                                  This should be a JSON file.
  --cap-size INTEGER              The maximum number of torsions to include
                                  per parameter in the auxiliary datasets.If
                                  there are more torsions than this, a subset
                                  will be selected.  [default: 5]
  --cap-method [pick_random|pick_heavy|pick_light]
                                  The method to use to select the torsions to
                                  include per parameter in the auxiliary
                                  datasets.  [default: pick_random]
  --verbose                       Whether to print out additional information.
  --n-processes INTEGER           The number of processes to use when
                                  processing the data.  [default: 4]
  --min-record-coverage INTEGER   The minimum number of records a parameter
                                  must have to be included in the force field
                                  optimization.  [default: 5]
  --help                          Show this message and exit.
```


