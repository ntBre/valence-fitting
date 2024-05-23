import tempfile

from openff.toolkit import ForceField, Molecule
from rdkit import Chem

from query import find_matches, into_params, load_want, mol_from_smiles
from store import (
    DBMol,
    Store,
    bits_to_elements,
    elements_to_bits,
    get_elements,
)


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


def test_load_want():
    got = len(load_want("want.params"))
    want = 51
    assert got == want


def test_elements_to_bits():
    mol = Chem.MolFromSmiles("CCO")
    got = elements_to_bits(get_elements(mol))
    want = 0b101000000
    assert got == want

    got = bits_to_elements(got)
    want = [6, 8]
    assert got == want


def test_store():
    with tempfile.NamedTemporaryFile() as f:
        s = Store(f.name)
        mol = DBMol("CC[At]", 3, 1 << 85 | 1 << 6)
        s.insert_molecule(mol)
        got = next(s.get_molecules()).get_elements()
        want = [6, 85]
        assert got == want


def test_query():
    from query import _main

    with tempfile.NamedTemporaryFile() as f:
        _main(8, 32, [], f.name, "openff-2.1.0.offxml", ["t1", "t2"], 100)
