# splitting more torsions with differing multiplicities. looking for smirks
# with commas and different XN values on either side. for example, #16X2,#16X1
# because those will cause different multiplicities. also only care about the
# central atoms since that's what determines the multiplicity

import re

from openff.toolkit import ForceField

# comma somewhere inside square brackets
pat = re.compile(r"\[([^]]+),([^]]+)\]")

ff = ForceField("output/initial-force-field-openff-2.1.0.offxml")
h = ff.get_parameter_handler("ProperTorsions")

for p in h.parameters:
    if s := pat.search(p.smirks):
        print(p.id, p.smirks)
