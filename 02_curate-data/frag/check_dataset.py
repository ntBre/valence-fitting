import logging
import re

from openff.toolkit import ForceField, Molecule

from serve import ffname

logging.getLogger("openff").setLevel(logging.ERROR)

ff = ForceField(ffname)

junk = re.compile("[(,)]")

conformers = 0
with open("dataset.smi") as inp:
    for line in inp:
        pid, cmiles, *rest = line.split()
        tors = tuple([int(junk.sub("", x)) for x in rest])
        mol = Molecule.from_mapped_smiles(cmiles, allow_undefined_stereo=True)
        print(pid, cmiles, tors)
        assert pid in [
            p.id
            for p in ff.label_molecules(mol.to_topology())[0][
                "ProperTorsions"
            ].values()
        ]
        try:
            mol.assign_partial_charges("am1bccelf10")
        except ValueError as e:
            if "Omega conformer generation failed" in str(e):
                conformers += 1
            else:
                raise e

print(f"{conformers} failed for omega errors")
