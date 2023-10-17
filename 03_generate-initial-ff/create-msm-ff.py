import faulthandler
import json
import logging
import os
import typing
from collections import defaultdict

import click
import numpy as np
import tqdm
from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResultCollection,
)
from openff.qcsubmit.results.filters import LowestEnergyFilter
from openff.toolkit import ForceField
from openff.units import unit
from qubekit.bonded.mod_seminario import ModSeminario, ModSemMaths
from qubekit.molecules import Ligand

if typing.TYPE_CHECKING:
    from openff.toolkit import Molecule
    from qcportal.models import ResultRecord

logging.getLogger("openff").setLevel(logging.ERROR)


# it's time for monkey patching
def force_constant_bond(bond, eigenvals, eigenvecs, coords):
    atom_a, atom_b = bond
    eigenvals_ab = eigenvals[atom_a, atom_b, :]
    eigenvecs_ab = eigenvecs[:, :, atom_a, atom_b]

    unit_vectors_ab = ModSemMaths.unit_vector_along_bond(coords, bond)

    return -0.5 * sum(
        eigenvals_ab[i] * abs(np.dot(unit_vectors_ab, eigenvecs_ab[:, i]))
        for i in range(3)
    )


ModSemMaths.force_constant_bond = force_constant_bond


def cross(u_ab, u_bc):
    "Return the square norm of the cross product of the arguments"
    c = np.cross(u_ab, u_bc)
    d = np.linalg.norm(c)
    return d * d


def calculate_torsions(eigenvals, eigenvecs, molecule, bond_lens):
    coords = molecule.coordinates
    ret = dict()
    # molecule.dihedrals is, for some unholy reason, a dict of central_atoms:
    # list_of_dihedrals_involving_them. so you have to call .values to get the
    # list and then iterate over the list to get actual dihedrals
    for dihedral in (
        dihedral
        for dihedrals in molecule.dihedrals.values()
        for dihedral in dihedrals
    ):
        atom_a, atom_b, atom_c, atom_d = dihedral
        u_ab = ModSemMaths.unit_vector_along_bond(coords, (atom_a, atom_b))
        u_cb = ModSemMaths.unit_vector_along_bond(coords, (atom_c, atom_b))
        u_bc = -u_cb
        u_dc = ModSemMaths.unit_vector_along_bond(coords, (atom_d, atom_b))
        u_cd = -u_dc

        u_n_abc = ModSemMaths.unit_vector_normal_to_bond(u_cb, u_ab)
        u_n_bcd = ModSemMaths.unit_vector_normal_to_bond(u_dc, u_bc)
        r_ba = bond_lens[atom_b, atom_a]
        eigenvals_ab = eigenvals[atom_a, atom_b, :]
        eigenvecs_ab = eigenvecs[:3, :3, atom_a, atom_b]

        r_cd = bond_lens[atom_c, atom_d]
        eigenvals_dc = eigenvals[atom_d, atom_c, :]
        eigenvecs_dc = eigenvecs[:3, :3, atom_d, atom_c]

        sum_first = sum(
            eigenvals_ab[i]
            * abs(ModSemMaths.dot_product(u_n_abc, eigenvecs_ab[:, i]))
            for i in range(3)
        )
        sum_second = sum(
            eigenvals_dc[i]
            * abs(ModSemMaths.dot_product(u_n_bcd, eigenvecs_dc[:, i]))
            for i in range(3)
        )
        k = 1.0 / (r_ba * r_ba * cross(u_ab, u_bc) * sum_first) + 1.0 / (
            r_cd * r_cd * cross(u_bc, u_cd) * sum_second
        )
        ret[dihedral] = np.real(1.0 / k)

    return ret


def _modified_seminario_method(
    self, molecule: Ligand, hessian: np.ndarray
) -> Ligand:
    """Calculate the new bond and angle terms after being passed the symmetric
    Hessian and optimised molecule coordinates.

    """
    size_mol = molecule.n_atoms
    eigenvecs = np.empty((3, 3, size_mol, size_mol), dtype=complex)
    eigenvals = np.empty((size_mol, size_mol, 3), dtype=complex)
    bond_lens = np.zeros((size_mol, size_mol))

    for i in range(size_mol):
        for j in range(size_mol):
            diff_i_j = molecule.coordinates[i, :] - molecule.coordinates[j, :]
            bond_lens[i, j] = np.linalg.norm(diff_i_j)

            partial_hessian = hessian[
                (i * 3) : ((i + 1) * 3), (j * 3) : ((j + 1) * 3)
            ]

            eigenvals[i, j, :], eigenvecs[:, :, i, j] = np.linalg.eig(
                partial_hessian
            )

    # The bond and angle values are calculated and written to file.
    self.calculate_bonds(eigenvals, eigenvecs, molecule, bond_lens)
    if molecule.n_angles > 0:
        # handle linear molecules with no angles
        self.calculate_angles(eigenvals, eigenvecs, molecule, bond_lens)
    torsions = []
    if molecule.n_dihedrals > 0:
        torsions = calculate_torsions(
            eigenvals, eigenvecs, molecule, bond_lens
        )
    return molecule, torsions


def run(self, molecule: Ligand, **kwargs) -> Ligand:
    import copy

    from qubekit.utils import constants

    # reset the bond and angle parameter groups
    molecule.BondForce.clear_parameters()
    molecule.AngleForce.clear_parameters()
    # convert the hessian from atomic units
    conversion = constants.HA_TO_KCAL_P_MOL / (constants.BOHR_TO_ANGS**2)
    # make sure we do not change the molecule hessian
    hessian = copy.deepcopy(molecule.hessian)
    hessian *= conversion
    # they don't even capture this return??
    _, torsions = self._modified_seminario_method(
        molecule=molecule, hessian=hessian
    )
    # apply symmetry to the bond and angle parameters
    molecule.symmetrise_bonded_parameters()

    return molecule, torsions


ModSeminario._modified_seminario_method = _modified_seminario_method
ModSeminario.run = run


def calculate_parameters(
    qc_record: "ResultRecord",
    molecule: "Molecule",
    forcefield: "ForceField",
) -> typing.Dict[str, typing.Dict[str, typing.List[unit.Quantity]]]:
    """
    Calculate the modified seminario parameters for the given input molecule
    and store them by OFF SMIRKS.
    """
    mod_sem = ModSeminario()

    # create the qube molecule, this should be in the same order as the off_mol
    qube_mol = Ligand.from_rdkit(molecule.to_rdkit(), name="offmol")
    qube_mol.hessian = qc_record.return_result
    # calculate the modified seminario parameters and store in the molecule
    # no idea what this is
    qube_mol, torsions = mod_sem.run(qube_mol)
    # label the openff molecule
    labels = forcefield.label_molecules(molecule.to_topology())[0]
    # loop over all bonds and angles and collect the results in
    # nm/ kj/mol / radians(openMM units)
    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
        "torsions": defaultdict(list),
    }

    for bond, parameter in labels["Bonds"].items():
        # bond is a tuple of the atom index the parameter is applied to
        qube_param = qube_mol.BondForce[bond]
        all_parameters["bond_eq"][parameter.smirks].append(qube_param.length)
        all_parameters["bond_k"][parameter.smirks].append(qube_param.k)

    for angle, parameter in labels["Angles"].items():
        qube_param = qube_mol.AngleForce[angle]
        all_parameters["angle_eq"][parameter.smirks].append(qube_param.angle)
        all_parameters["angle_k"][parameter.smirks].append(qube_param.k)

    for torsion, parameter in labels["ProperTorsions"].items():
        qube_param = torsions.get(torsion, None)
        if qube_param is not None:
            all_parameters["torsions"][parameter.smirks].append(qube_param)

    return all_parameters


@click.command()
@click.option(
    "--initial-force-field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the initial force field file (OFFXML).",
)
@click.option(
    "--output",
    "output_force_field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file (OFFXML).",
)
@click.option(
    "--optimization-dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset.",
)
@click.option(
    "--working-directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=False,
    help=(
        "The path to the working directory. "
        "Intermediate files are saved here if provided"
    ),
)
@click.option(
    "--verbose/--no-verbose",
    default=False,
    help="Enable verbose logging.",
)
def main(
    initial_force_field: str,
    output_force_field: str,
    optimization_dataset: str,
    working_directory: typing.Optional[str] = None,
    verbose: bool = False,
):
    dataset = OptimizationResultCollection.parse_file(optimization_dataset)

    # filter for lowest energy results
    print("filtering")
    filtered = dataset.filter(LowestEnergyFilter())

    # filter to only keep entries with hessians calculated
    print("converting to results")
    hessian_set = filtered.to_basic_result_collection(driver="hessian")

    if working_directory is not None:
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)
        hessian_file = os.path.join(working_directory, "hessian_set.json")
        with open(hessian_file, "w") as f:
            f.write(hessian_set.json(indent=2))
        if verbose:
            print(f"Hessian set written to: {hessian_file}")

    if verbose:
        print(f"Found {hessian_set.n_results} hessian calculations")
        print(f"Found {hessian_set.n_molecules} hessian molecules")

    ff = ForceField(initial_force_field, allow_cosmetic_attributes=True)

    records_and_molecules = list(hessian_set.to_records())
    if verbose:
        records_and_molecules = tqdm.tqdm(
            records_and_molecules,
            desc="Calculating parameters",
        )

    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
        "torsions": defaultdict(list),
    }
    errored_records_and_molecules = []
    for record, molecule in records_and_molecules:
        try:
            parameters = calculate_parameters(record, molecule, ff)
        except BaseException:
            errored_records_and_molecules.append((record, molecule))
            continue
        else:
            for key, values in parameters.items():
                for smirks, value in values.items():
                    all_parameters[key][smirks].extend(value)

    if working_directory is not None:
        seminario_file = os.path.join(
            working_directory, "seminario_parameters.json"
        )
        with open(seminario_file, "w") as file:
            json.dump(all_parameters, file, indent=2)

    if verbose:
        print(
            f"Found {len(errored_records_and_molecules)} errored calculations"
        )
    if working_directory is not None:
        if len(errored_records_and_molecules):
            key = list(dataset.entries.keys())[0]
            opt_records_by_id = {
                record.record_id: record for record in hessian_set.entries[key]
            }
            records, _ = zip(*errored_records_and_molecules)
            errored_records = [
                opt_records_by_id[record.id] for record in records
            ]
            errored_dataset = BasicResultCollection(
                entries={key: errored_records}
            )
            error_file = os.path.join(
                working_directory, "errored_dataset.json"
            )
            with open(error_file, "w") as f:
                f.write(errored_dataset.json(indent=2))
            if verbose:
                print(f"Errored dataset written to: {error_file}")

    # now we need to update the FF parameters
    kj_per_mol_per_nm2 = unit.kilojoule_per_mole / unit.nanometer**2
    bond_handler = ff.get_parameter_handler("Bonds")
    for smirks in all_parameters["bond_eq"]:
        bond = bond_handler.parameters[smirks]

        bond_length = (
            np.mean(all_parameters["bond_eq"][smirks]) * unit.nanometer
        )
        bond.length = bond_length.to(unit.angstrom)

        bond_k = np.mean(all_parameters["bond_k"][smirks]) * kj_per_mol_per_nm2
        bond.k = bond_k.to(unit.kilocalorie_per_mole / (unit.angstrom**2))

    kj_per_mol_per_rad2 = unit.kilojoule_per_mole / (unit.radian**2)
    angle_handler = ff.get_parameter_handler("Angles")
    for smirks in all_parameters["angle_eq"]:
        angle = angle_handler.parameters[smirks]

        angle_eq = np.mean(all_parameters["angle_eq"][smirks]) * unit.radian
        angle.angle = angle_eq.to(unit.degree)

        angle_k = (
            np.mean(all_parameters["angle_k"][smirks]) * kj_per_mol_per_rad2
        )
        angle.k = angle_k.to(unit.kilocalorie_per_mole / unit.radian**2)

    torsion_handler = ff.get_parameter_handler("ProperTorsions")
    for smirks, vals in all_parameters["torsions"].items():
        # TODO figure out these units. I'm assuming kJ like the others and the
        # denominator has four distance terms in it, which should give nm^4,
        # but our force constants don't have distance units
        val = np.mean(vals) * unit.kilojoule_per_mole  # / unit.nanometer**4
        print(torsion_handler[smirks].k1, val.to(unit.kilocalorie_per_mole))

    ff.to_file(output_force_field)


if __name__ == "__main__":
    with open("fault_handler.log", "w") as fobj:
        faulthandler.enable(fobj)
        main()
