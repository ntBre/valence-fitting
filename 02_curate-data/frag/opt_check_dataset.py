import logging

from openff.toolkit import ForceField, Molecule

from serve import ffname

logging.getLogger("openff").setLevel(logging.ERROR)

ff = ForceField(ffname)


conformers = []
with open("remaining.opt") as inp:
    for line in inp:
        if line.startswith("#"):
            continue
        pid, _smarts, cmiles = line.split()
        mol = Molecule.from_smiles(cmiles, allow_undefined_stereo=True)
        print(pid, cmiles)

        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        pids = [p.id for p in labels.values()]
        assert pid in pids
        try:
            mol.assign_partial_charges("am1bccelf10")
        except ValueError as e:
            if "Omega conformer generation failed" in str(e):
                conformers.append(cmiles)
            else:
                raise e

print(f"{len(conformers)} failed for omega errors")
for conf in conformers:
    print(conf)
