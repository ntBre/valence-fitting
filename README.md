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