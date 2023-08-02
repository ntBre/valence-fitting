# Generating initial force field

In this directory we generate the initial force field to re-fit. If parameter
SMIRKS need to be modified, that should happen here. This will provide the input
to steps 2 (curating data for fitting these SMIRKS) and 3 (generating initial
parameter values).

In this case, `generate_forcefield.py` computes the differences between the
original Sage 2.0 force field and the updated version with additional torsions
(`force-field.offxml`), and ports these changes to Sage 2.1. The resulting force
field is Sage 2.1 with the additional torsions.
