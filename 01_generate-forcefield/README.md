# Generating initial force field

In this directory we generate the initial force field to re-fit. If parameter SMIRKS need to be modified, that should happen here. This will provide the input to steps 2 (curating data for fitting these SMIRKS) and 3 (generating initial parameter values).

The `generate_forcefield.py` script here is just a placeholder that should be replaced if the parameters are modified. It downloads a previous force field for re-fitting.

```
$ ./generate-forcefield.sh
```