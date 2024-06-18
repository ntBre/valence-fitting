import logging
import re
import sys

from openff.toolkit import ForceField, Molecule

from serve import ffname

logging.getLogger("openff").setLevel(logging.ERROR)

ff = ForceField(ffname)

junk = re.compile("[(,)]")

conformers = []
with open(sys.argv[1]) as inp:
    for line in inp:
        if line.startswith("#"):
            continue
        pid, cmiles, *rest = line.split()
        tors = tuple([int(junk.sub("", x)) for x in rest])
        mol = Molecule.from_mapped_smiles(cmiles, allow_undefined_stereo=True)
        print(pid, cmiles, tors)

        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        assert labels[tors].id == pid
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
