# Valence fitting

This repository is a small modification of scripts at
https://github.com/openforcefield/sage-2.1.0.

Create the environment necessary with the below. Replace ``micromamba`` with
``mamba`` or ``conda`` as needed.

```bash
micromamba env create -f environment.yaml
```

If you're running on Apple Silicon, some packages are still not available for
ARM yet. Instead, do the following:


```bash
CONDA_SUBDIR=osx-64 micromamba env create -f environment.yaml
micromamba activate fb-195-tk-014-py310
conda config --env --set subdir osx-64  # make sure updates keep using the subdir
```

The general steps laid out here are:

1. Generate the initial force field for fitting. This step defines the SMIRKS of
   each parameter that will be fit. When adding, removing, or modifying
   parameters in force fields, here is where it should be changed.
2. Download and curate QM datasets. This step downloads QM datasets; manually
   adds additional data records; and manually removes bad records. For torsions,
   you can specify whether to cap the number of torsion data records per torsion
   parameter for fitting. For optimizations, you can cap the number of
   conformers used as well. The output of these scripts are the QM datasets, as
   well as which parameters will be trained. A parameter will not be trained if
   there isn't enough data coverage (also specified by the user).
3. This step generates initial values for each force field parameter using the
   Modified Seminario Method.
4. This step sets up the final ForceBalance input files for training.

After setting up ForceBalance, the ``scripts`` directory contains scripts
suitable for running the fit on HPC3.


The included `Makefile` can be used to run the whole pipeline at once:

```bash
make step4
```

Alternatively, one of the preceding steps can be run directly. `make` will
handle generating any required dependencies. For example,

```bash
make step3
```

will generate `03_generate-initial-ff/output/initial-force-field-msm.offxml`,
which depends on
`01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml` from
step 1, `02_curate-data/output/combined-opt.json` from step 2, and the script
`03_generate-initial-ff/create-msm-ff.py` from step 3. These files from steps 1
and 2 in turn depend on other input files and scripts, which will be rerun if
their own dependencies change.
