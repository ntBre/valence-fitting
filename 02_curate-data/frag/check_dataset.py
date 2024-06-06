import re

from openff.toolkit import ForceField, Molecule

from serve import ffname

ff = ForceField(ffname)

junk = re.compile("[(,)]")

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
