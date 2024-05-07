from openff.toolkit import ForceField, Molecule

from query import find_matches, into_params, mol_from_smiles


def test_find_matches():
    ff = ForceField("openff-2.1.0.offxml")
    params = into_params(ff)

    for smile in ["CCO"]:
        offmol = Molecule.from_smiles(smile)
        rdmol = mol_from_smiles(smile)
        got = list(sorted(find_matches(params, rdmol).items()))
        want = list(
            (k, v.id)
            for k, v in sorted(
                ff.label_molecules(offmol.to_topology())[0][
                    "ProperTorsions"
                ].items()
            )
        )
        assert got == want, f"got = {got}, want = {want}"
