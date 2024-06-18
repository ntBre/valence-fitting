from openff.toolkit import ForceField

from serve import ffname

with open("want.opt") as inp:
    want = {line.strip() for line in inp}

with open("opt.smi") as inp:
    got = {line.strip().split()[0] for line in inp}


ff = ForceField(ffname)
tors = ff.get_parameter_handler("ProperTorsions")
with open("remaining.opt", "w") as out:
    for p in want - got:
        print(p, tors.get_parameter(dict(id=p))[0].smirks, file=out)
