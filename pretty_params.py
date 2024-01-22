from openff.toolkit import ForceField

with open("interesting_params.txt") as inp:
    params = ["t" + line.strip() for line in inp]

h = ForceField("openff-2.1.0.offxml").get_parameter_handler("ProperTorsions")

for pid in params:
    print(f"{pid:>4s} {h.get_parameter(dict(id=pid))[0].smirks}")
