from collections import defaultdict
from copy import deepcopy

from openff.toolkit import ForceField

base = ForceField("openff-2.1.0.offxml")

with open("to_add.dat") as inp:
    to_add = [
        line.strip()
        for line in inp
        if not line.startswith("#") and len(line) > 1
    ]

new_params = defaultdict(list)
for line in to_add:
    [pid, smirks] = line.split()
    new_params[pid].append(smirks)

to_remove = ["t123"]

torsions = base.get_parameter_handler("ProperTorsions")

old = deepcopy(torsions.parameters)
torsions.parameters.clear()

for t in old:
    if t.id in new_params:
        for i, p in enumerate(new_params[t.id]):
            tmp = deepcopy(t)
            tmp.id = t.id + "ghijkl"[i]
            tmp.smirks = p
            torsions.parameters.append(tmp)
    elif t.id not in to_remove:
        torsions.parameters.append(t)

base.to_file("output/initial-force-field-openff-2.1.0.offxml")
