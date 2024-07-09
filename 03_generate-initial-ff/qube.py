import abc
import copy
import decimal
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from enum import Enum
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple, Union
from xml.dom.minidom import parseString

import networkx as nx
import numpy as np
import qcelemental as qcel
import qcengine as qcng
import qubekit
from openff.toolkit.typing.engines.smirnoff import get_available_force_fields
from pydantic.v1 import BaseModel, Field, PositiveInt, dataclasses, validator
from qcelemental.models.types import Array
from qubekit.forcefield import (
    BaseForceGroup,
    HarmonicAngleForce,
    HarmonicBondForce,
    LennardJones126Force,
    PeriodicImproperTorsionForce,
    PeriodicTorsionForce,
    RBImproperTorsionForce,
    RBProperTorsionForce,
    VirtualSiteGroup,
)
from qubekit.utils.exceptions import (
    ConformerError,
    FileTypeError,
    MissingReferenceData,
    SmartsError,
    SpecificationError,
    StereoChemistryError,
    TopologyMismatch,
    TorsionDriveDataError,
)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.rdchem import GetPeriodicTable, PeriodicTable
from rdkit.Geometry.rdGeometry import Point3D
from typing_extensions import Literal

if TYPE_CHECKING:
    from qubekit.molecules import Ligand

try:
    # fix for the openmm namechange
    from openmm import unit
    from openmm.app import Aromatic, Double, PDBFile, Single, Topology, Triple
    from openmm.app.element import Element
except (ModuleNotFoundError, ImportError):
    from simtk import unit
    from simtk.openmm.app import (
        Aromatic,
        Double,
        PDBFile,
        Single,
        Topology,
        Triple,
    )
    from simtk.openmm.app.element import Element

AVOGADRO = 6.02214179e23  # Particles in 1 Mole
ROOM_TEMP = 298.15  # Kelvin
ROOM_PRESSURE = 101325  # Pascals
VACUUM_PERMITTIVITY = 8.8541878128e-12  # Farads per Metre
ELECTRON_CHARGE = 1.602176634e-19  # Coulombs
KB_KCAL_P_MOL_K = 0.0019872041  # Boltzmann constant in KCal/(mol * K)

PI = 3.141592653589793  # Pi
DEG_TO_RAD = PI / 180  # Degrees to radians
RAD_TO_DEG = 180 / PI  # Radians to degrees


KCAL_TO_KJ = 4.184  # Kilocalories to kiloJoules
KJ_TO_KCAL = 0.23900573613  # KiloJoules to kilocalories
J_TO_KCAL = 0.0002390057  # Joules to kilocalories

J_TO_KCAL_P_MOL = J_TO_KCAL * AVOGADRO  # Joules to kilocalories per mole

HA_TO_KCAL_P_MOL = 627.509391  # Hartrees to kilocalories per mole
KCAL_P_MOL_TO_HA = 0.00159360164  # Kilocalories per mole to Hartrees

NM_TO_ANGS = 10  # Nanometres to Angstroms
ANGS_TO_NM = 0.1  # Angstroms to nanometres

ANGS_TO_M = 1e-10  # Angstroms to metres
M_TO_ANGS = 1e10  # Metres to Angstroms

BOHR_TO_ANGS = 0.529177  # Bohrs to Angstroms
ANGS_TO_BOHR = 1.88972687777  # Angstroms to Bohrs

EPSILON_CONVERSION = (
    (BOHR_TO_ANGS**6) * HA_TO_KCAL_P_MOL * KCAL_TO_KJ
)  # L-J Conversion
SIGMA_CONVERSION = ANGS_TO_NM  # L-J Conversion


class SchemaBase(BaseModel):
    """A basic pydantic starting class which uses assigment validation."""

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True
        json_encoders = {np.ndarray: lambda v: v.flatten().tolist()}


class TDSettings(SchemaBase):
    """
    A schema with available options for Time-Dependent calculations.
    """

    n_states: int = Field(3, description="The number of states to solve for.")
    use_tda: bool = Field(
        False,
        description="If we should use the Tamm-Dancoff approximation (TDA).",
    )


class QCOptions(SchemaBase):
    """
    A simple Schema to validate QC/ML/MM runtime options.
    Note this model is locked once created to avoid validation errors.
    """

    program: str = Field(
        "gaussian",
        description="The name of the program which should be used to carry out the computation, such as psi4",
    )
    basis: Optional[str] = Field(
        "6-311++G(d,p)",
        description="The basis that should be used in the computation.",
    )
    method: str = Field(
        "wB97X-D",
        description="The method that should be used for the computation.",
    )
    td_settings: Optional[TDSettings] = Field(
        None,
        description="Any time dependent settings that should be used during the computation. Note not all programs support this option.",
    )

    @validator("program", "method")
    def _cast_lower(cls, parameter: str) -> str:
        """Lower the parameter to avoid validation issues."""
        return parameter.lower()

    def validate_program(self):
        """
        Validate the choice of program against those supported by QCEngine and QUBEKit.
        """
        programs = qcng.list_available_programs()
        programs.discard("dftd3")

        if self.program.lower() not in programs:
            raise SpecificationError(
                f"The program {self.program} is not available, available programs are {programs}"
            )

    @property
    def keywords(self) -> Dict[str, Union[str, int]]:
        """
        Build some keywords in a consistent way for the qcspec.
        """
        keywords = {
            "scf_type": "df",
            # make sure we always use an ultrafine grid
            "dft_spherical_points": 590,
            "dft_radial_points": 99,
        }
        if self.td_settings is not None:
            # use psi4 keyword settings to be consistent
            keywords["tdscf_states"] = self.td_settings.n_states
            keywords["tdscf_tda"] = self.td_settings.use_tda

        # work around a setting in psi4, fixes range seperated functionals
        if self.program.lower() == "psi4":
            keywords["wcombine"] = False
        return keywords

    @property
    def qc_model(self) -> qcel.models.common_models.Model:
        """
        Build the QC model for the computation.

        Important:
            The method name can be changed depending on the program used and td settings
        """
        if self.td_settings is not None and self.program == "psi4":
            # we have to add the td tag
            method = self.method
            if "td" != method.split("-")[0]:
                method = f"td-{method}"
        else:
            method = self.method
        model = qcel.models.common_models.Model(
            method=method, basis=self.basis
        )
        return model

    def validate_specification(self) -> None:
        """
        Validate the specification this should be called before using the spec to find errors.
        """
        # make sure the program is valid first then the basis method combination
        self.validate_program()

        openff_forcefields = [
            ff.split(".offxml")[0].lower()
            for ff in get_available_force_fields()
        ]
        # set up some models
        ani_methods = {"ani1x", "ani1ccx", "ani2x"}
        xtb_methods = {
            "gfn0-xtb",
            "gfn0xtb",
            "gfn1-xtb",
            "gfn1xtb",
            "gfn2-xtb",
            "gfn2xtb",
            "gfn-ff",
            "gfnff",
        }
        rdkit_methods = {"uff", "mmff94", "mmff94s"}
        gaff_forcefields = {
            "gaff-1.4",
            "gaff-1.8",
            "gaff-1.81",
            "gaff-2.1",
            "gaff-2.11",
        }
        settings = {
            "openmm": {
                "antechamber": gaff_forcefields,
                "smirnoff": openff_forcefields,
            },
            "torchani": {None: ani_methods},
            "xtb": {None: xtb_methods},
            "rdkit": {None: rdkit_methods},
        }
        # now check these settings
        # TODO do we raise an error or just change at run time with a warning?
        # we do not validate QM as there are so many options
        if self.program.lower() in settings:
            program_settings = settings[self.program.lower()]

            allowed_methods = program_settings.get(self.basis, None)
            if allowed_methods is None:
                raise SpecificationError(
                    f"The Basis {self.basis} is not supported for the program {self.program} please chose from {program_settings.keys()}"
                )
            # now check the method
            method = self.method.split(".offxml")[0].lower()
            if method not in allowed_methods:
                raise SpecificationError(
                    f"The method {method} is not available for the program {self.program}  with basis {self.basis}, please chose from {allowed_methods}"
                )
        if self.td_settings is not None:
            if self.program.lower() not in ["gaussian"]:
                raise SpecificationError(
                    f"The program {self.program.lower()} does not support time-dependent calculations."
                )


class StageBase(SchemaBase, abc.ABC):

    type: Literal["base"] = "base"

    @classmethod
    @abc.abstractmethod
    def is_available(cls) -> bool:
        """Check any dependencies to make sure that this stage is available to run."""
        ...

    @abc.abstractmethod
    def run(self, molecule: "Ligand", **kwargs) -> "Ligand":
        """The main function of the stage which should perform some parametrisation and return the complete molecule."""
        ...

    @abc.abstractmethod
    def start_message(self, **kwargs) -> str:
        """
        A friendly message to let users know that stage is starting with any important options.
        """
        ...

    @abc.abstractmethod
    def finish_message(self, **kwargs) -> str:
        """
        A friendly message to let users know that the stage is complete and any checks that have been performed.
        """
        ...


@dataclasses.dataclass
class TorsionScan:
    torsion: Tuple[int, int, int, int]
    scan_range: Tuple[int, int]


@dataclasses.dataclass
class GridPointResult:
    """A simple class to help with the torsiondrive json API.
    Important:
        geometries are in bohr
        energy in hartree
    """

    dihedral_angle: int
    input_geometry: List[float]
    final_geometry: List[float]
    final_energy: float


class AtomStereoChemistry(str, Enum):
    """
    Atom stereochemistry types.
    """

    R = "R"
    S = "S"
    U = "Unknown"


class BondStereoChemistry(str, Enum):
    """
    Bond stereochemistry types.
    """

    E = "E"
    Z = "Z"
    U = "Unknown"


@dataclasses.dataclass  # Cannot be frozen as params are loaded separately.
class AIM:
    """Data in atomic units."""

    volume: Optional[float] = None
    charge: Optional[float] = None
    c8: Optional[float] = None
    # TODO Extend to include other types of potential e.g. Buckingham


@dataclasses.dataclass(frozen=True)
class Dipole:
    """Data in atomic units."""

    x: float
    y: float
    z: float

    def to_array(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclasses.dataclass(frozen=True)
class Quadrupole:
    """Data in atomic units."""

    q_xx: float
    q_xy: float
    q_xz: float
    q_yz: float
    q_yy: float
    q_zz: float

    def to_array(self) -> np.ndarray:
        return np.array(
            [
                [self.q_xx, self.q_xy, self.q_xz],
                [self.q_xy, self.q_yy, self.q_yz],
                [self.q_xz, self.q_yz, self.q_zz],
            ]
        )


@dataclasses.dataclass(frozen=True)
class CloudPen:
    """Data in atomic units."""

    a: float
    b: float


class Atom(BaseModel):
    """
    Class to hold all of the "per atom" information.
    All atoms in Molecule will have an instance of this Atom class to describe their properties.
    """

    class Config:
        validate_assignment = True
        json_encoders = {Enum: lambda v: v.value}

    atomic_number: int = Field(
        ...,
        description="The atomic number of the atom all other properties are based on this number.",
        ge=0,
    )
    atom_index: int = Field(
        ...,
        description="The index this atom has in the molecule object",
        ge=0,
    )
    atom_name: Optional[str] = Field(
        None,
        description="An optional unqiue atom name that should be assigned to the atom, the ligand object will make sure all atoms have unique names.",
    )
    formal_charge: int = Field(
        ...,
        description="The formal charge of the atom, used to calculate the molecule total charge",
    )
    aromatic: bool = Field(
        ...,
        description="If the atom should be considered aromatic `True` or not `False`.",
    )
    stereochemistry: Optional[AtomStereoChemistry] = Field(
        None,
        description="The stereochemistry of the atom where None means not stereogenic and U is unknown or ambiguous.",
    )
    bonds: Optional[List[int]] = Field(
        None,
        description="The list of atom indices which are bonded to this atom.",
    )
    aim: Optional[AIM] = Field(
        AIM(),
    )
    dipole: Optional[Dipole] = Field(
        None,
    )
    quadrupole: Optional[Quadrupole] = Field(
        None,
    )
    cloud_pen: Optional[CloudPen] = Field(
        None,
    )

    @classmethod
    def from_rdkit(cls, rd_atom: Chem.Atom) -> "Atom":
        """
        Build a QUBEKit atom from an rdkit atom instance.
        """
        atomic_number = rd_atom.GetAtomicNum()
        index = rd_atom.GetIdx()
        formal_charge = rd_atom.GetFormalCharge()
        aromatic = rd_atom.GetIsAromatic()
        bonds = [a.GetIdx() for a in rd_atom.GetNeighbors()]
        # check for names in the normal places pdb, mol2 and mol
        if rd_atom.HasProp("_Name"):
            name = rd_atom.GetProp("_Name")
        elif rd_atom.HasProp("_TriposAtomName"):
            name = rd_atom.GetProp("_TriposAtomName")
        else:
            try:
                name = rd_atom.GetMonomerInfo().GetName().strip()
            except AttributeError:
                name = None
        # stereochem
        if rd_atom.HasProp("_CIPCode"):
            stereo_code = rd_atom.GetProp("_CIPCode")
        else:
            stereo_code = None
        return cls(
            atomic_number=atomic_number,
            atom_index=index,
            atom_name=name,
            formal_charge=formal_charge,
            aromatic=aromatic,
            stereochemistry=stereo_code,
            bonds=bonds,
        )

    @property
    def atomic_mass(self) -> float:
        """Convert the atomic number to mass."""
        return Element.mass(self.atomic_number)

    @property
    def atomic_symbol(self) -> str:
        """Convert the atomic number to the atomic symbol as per the periodic table."""
        return Element.name(self.atomic_number).title()

    def to_rdkit(self) -> Chem.Atom:
        """
        Convert the QUBEKit atom an RDKit atom.
        """
        # build the atom from atomic number
        rd_atom = Chem.Atom(self.atomic_number)
        rd_atom.SetFormalCharge(self.formal_charge)
        rd_atom.SetIsAromatic(self.aromatic)
        rd_atom.SetProp("_Name", self.atom_name)
        # left is counter clockwise
        if self.stereochemistry == "S":
            rd_atom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
        # right is clockwise
        elif self.stereochemistry == "R":
            rd_atom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)

        return rd_atom


class Bond(BaseModel):
    """
    A basic bond class.
    """

    class Config:
        validate_assignment = True
        json_encoders = {Enum: lambda v: v.value}

    atom1_index: int = Field(
        ..., description="The index of the first atom in the bond."
    )
    atom2_index: int = Field(
        ..., description="The index of the second atom in the bond."
    )
    bond_order: float = Field(
        ..., description="The float value of the bond order."
    )
    aromatic: bool = Field(
        ..., description="If the bond should be considered aromatic."
    )
    stereochemistry: Optional[BondStereoChemistry] = Field(
        None,
        description="The stereochemistry of the bond, where None means not stereogenic.",
    )

    @classmethod
    def from_rdkit(cls, rd_bond: Chem.Bond) -> "Bond":
        """
        Build a QUBEKit bond class from an rdkit reference.
        """
        atom1_index = rd_bond.GetBeginAtomIdx()
        atom2_index = rd_bond.GetEndAtomIdx()
        aromatic = rd_bond.GetIsAromatic()
        order = rd_bond.GetBondTypeAsDouble()
        stereo_tag = rd_bond.GetStereo()
        if stereo_tag == Chem.BondStereo.STEREOZ:
            stereo = "Z"
        elif stereo_tag == Chem.BondStereo.STEREOE:
            stereo = "E"
        else:
            stereo = None
        return cls(
            atom1_index=atom1_index,
            atom2_index=atom2_index,
            aromatic=aromatic,
            bond_order=order,
            stereochemistry=stereo,
        )

    @property
    def rdkit_type(self) -> Chem.BondType:
        """
        Convert the bond order float to a bond type.
        """
        conversion = {
            1: Chem.BondType.SINGLE,
            1.5: Chem.BondType.AROMATIC,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.QUADRUPLE,
            5: Chem.BondType.QUINTUPLE,
            6: Chem.BondType.HEXTUPLE,
            7: Chem.BondType.ONEANDAHALF,
        }
        return conversion[self.bond_order]

    @property
    def rdkit_stereo(self) -> Optional[Chem.BondStereo]:
        """
        Return the rdkit style stereo enum.
        """
        if self.stereochemistry == "E":
            return Chem.BondStereo.STEREOE
        elif self.stereochemistry == "Z":
            return Chem.BondStereo.STEREOZ
        return None

    @property
    def indices(self) -> Tuple[int, int]:
        return self.atom1_index, self.atom2_index


class RDKit:
    """Class for controlling useful RDKit functions."""

    @staticmethod
    def mol_to_file(rdkit_mol: Chem.Mol, file_name: str) -> None:
        """
        Write the rdkit molecule to the requested file type.
        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
            file_name:
                Name of the file to be created.
        """
        file_path = Path(file_name)
        if file_path.suffix == ".pdb":
            return Chem.MolToPDBFile(rdkit_mol, file_name)
        elif file_path.suffix == ".sdf" or file_path.suffix == ".mol":
            return Chem.MolToMolFile(rdkit_mol, file_name)
        elif file_path.suffix == ".xyz":
            return Chem.MolToXYZFile(rdkit_mol, file_name)
        else:
            raise FileTypeError(
                f"The file type {file_path.suffix} is not supported please chose from xyz, pdb, mol or sdf."
            )

    @staticmethod
    def mol_to_multiconformer_file(
        rdkit_mol: Chem.Mol, file_name: str
    ) -> None:
        """
        Write the rdkit molecule to a multi conformer file.
        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
            file_name:
                Name of the file to be created.
        """
        file_path = Path(file_name)
        # get the file block writer
        if file_path.suffix == ".pdb":
            writer = Chem.MolToPDBBlock
        elif file_path.suffix == ".mol" or file_path.suffix == ".sdf":
            writer = Chem.MolToMolBlock
        elif file_path.suffix == ".xyz":
            writer = Chem.MolToXYZBlock
        else:
            raise FileTypeError(
                f"The file type {file_path.suffix} is not supported please chose from xyz, pdb, mol or sdf."
            )
        with open(file_name, "w") as out:
            for i in range(rdkit_mol.GetNumConformers()):
                out.write(writer(rdkit_mol, confId=i))

    @staticmethod
    def file_to_rdkit_mol(file_path: Path) -> Chem.Mol:
        """
        Args:
            file_path:
                Path of the file used to generate the rdkit molecule.
        return:
            RDKit molecule object generated from its file (or None if incorrect file type is provided).
        """

        # Read the file
        if file_path.suffix == ".pdb":
            mol = Chem.MolFromPDBFile(
                file_path.as_posix(), removeHs=False, sanitize=False
            )
        elif file_path.suffix == ".mol2":
            mol = Chem.MolFromMol2File(
                file_path.as_posix(), removeHs=False, sanitize=False
            )
        elif file_path.suffix == ".mol" or file_path.suffix == ".sdf":
            mol = Chem.MolFromMolFile(
                file_path.as_posix(),
                removeHs=False,
                sanitize=False,
                strictParsing=True,
            )
        else:
            raise FileTypeError(
                f"The file type {file_path.suffix} is not supported."
            )
        # run some sanitation
        Chem.SanitizeMol(
            mol,
            (
                Chem.SANITIZE_ALL
                ^ Chem.SANITIZE_SETAROMATICITY
                ^ Chem.SANITIZE_ADJUSTHS
            ),
        )
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        Chem.AssignStereochemistryFrom3D(mol)

        # set the name of the input file
        mol.SetProp("_Name", file_path.stem)

        return mol

    @staticmethod
    def smiles_to_rdkit_mol(smiles_string: str, name: Optional[str] = None):
        """
        Converts smiles strings to RDKit mol object.
        Args:
            smiles_string:
                The hydrogen free smiles string
            name:
                The name of the molecule this will be used when writing the pdb file
        return:
            The RDKit molecule
        """

        mol = AllChem.MolFromSmiles(smiles_string)
        if name is None:
            name = input("Please enter a name for the molecule:\n>")
        mol.SetProp("_Name", name)
        mol_hydrogens = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol_hydrogens, randomSeed=1)
        AllChem.SanitizeMol(mol_hydrogens)
        return mol_hydrogens

    @staticmethod
    def rdkit_descriptors(rdkit_mol: Chem.Mol) -> Dict[str, float]:
        """
        Use RDKit Descriptors to extract properties and store in Descriptors dictionary.
        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
        returns:
            descriptors dictionary
        """

        # Use RDKit Descriptors to extract properties and store in Descriptors dictionary
        return {
            "Heavy atoms": Descriptors.HeavyAtomCount(rdkit_mol),
            "H-bond donors": Descriptors.NumHDonors(rdkit_mol),
            "H-bond acceptors": Descriptors.NumHAcceptors(rdkit_mol),
            "Molecular weight": Descriptors.MolWt(rdkit_mol),
            "LogP": Descriptors.MolLogP(rdkit_mol),
        }

    @staticmethod
    def get_smiles(
        rdkit_mol: Chem.Mol,
        isomeric: bool = True,
        explicit_hydrogens: bool = True,
        mapped: bool = False,
    ) -> str:
        """
        Use RDKit to generate a smiles string for the molecule.

        We work with a copy of the input molecule as we may assign an atom map number which
        will affect the CIP algorithm and could break symmetry groups.

        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
            isomeric:
                If True, the smiles should encode the stereochemistry.
            explicit_hydrogens:
                If True, the smiles should explicitly encode hydrogens.
            mapped:
                If True, the smiles should be mapped to preserve the ordering of the molecule.

        Returns:
            A string which encodes the molecule smiles corresponding the the input options.
        """
        cp_mol = copy.deepcopy(rdkit_mol)
        if mapped:
            explicit_hydrogens = True
            for atom in cp_mol.GetAtoms():
                # mapping starts from 1 as 0 means no mapping in rdkit
                atom.SetAtomMapNum(atom.GetIdx() + 1)
        if not explicit_hydrogens:
            cp_mol = Chem.RemoveHs(cp_mol)
        return Chem.MolToSmiles(
            cp_mol, isomericSmiles=isomeric, allHsExplicit=explicit_hydrogens
        )

    @staticmethod
    def get_smirks_matches(
        rdkit_mol: Chem.Mol, smirks: str
    ) -> List[Tuple[int, ...]]:
        """
        Query the molecule for the tagged smarts pattern (OpenFF SMIRKS).

        Args:
            rdkit_mol:
                The rdkit molecule instance that should be checked against the smarts pattern.
            smirks:
                The tagged SMARTS pattern that should be checked against the molecule.

        Returns:
            A list of atom index tuples which match the corresponding tagged atoms in the smarts pattern.
            Note only tagged atoms indices are returned.
        """
        cp_mol = copy.deepcopy(rdkit_mol)
        smarts_mol = Chem.MolFromSmarts(smirks)
        if smarts_mol is None:
            raise SmartsError(
                f"RDKit could not understand the query {smirks} please check again."
            )
        # we need a mapping between atom map and index in the smarts mol
        # to work out the index of the matched atom
        mapping = {}
        for atom in smarts_mol.GetAtoms():
            smart_index = atom.GetAtomMapNum()
            if smart_index != 0:
                # atom was tagged in the smirks
                mapping[smart_index - 1] = atom.GetIdx()
        # smarts can match forward and backwards so condense the matches
        all_matches = set()
        for match in cp_mol.GetSubstructMatches(
            smarts_mol, uniquify=True, useChirality=True
        ):
            smirks_atoms = [match[atom] for atom in mapping.values()]
            all_matches.add(tuple(smirks_atoms))
        return list(all_matches)

    @staticmethod
    def get_smarts(rdkit_mol: Chem.Mol) -> str:
        """
        Use RDKit to get the smarts string of the molecule.
        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
        return:
            The smarts string of the molecule
        """

        return Chem.MolToSmarts(rdkit_mol)

    @staticmethod
    def generate_conformers(
        rdkit_mol: Chem.Mol, conformer_no: int
    ) -> List[np.ndarray]:
        """
        Generate a set of conformers for the molecule including the input conformer.
        Args:
            rdkit_mol:
                A complete Chem.Mol instance of a molecule.
            conformer_no:
                The number of conformers made for the molecule
        return:
            A list of conformer position arrays
        """

        AllChem.EmbedMultipleConfs(
            rdkit_mol,
            numConfs=conformer_no,
            randomSeed=1,
            clearConfs=False,
            useBasicKnowledge=True,
            pruneRmsThresh=1,
            enforceChirality=True,
        )
        positions = rdkit_mol.GetConformers()

        return [conformer.GetPositions() for conformer in positions]

    @staticmethod
    def find_symmetry_classes(rdkit_mol: Chem.Mol) -> Dict[int, str]:
        """
        Generate list of tuples of symmetry-equivalent (homotopic) atoms in the molecular graph
        based on: https://sourceforge.net/p/rdkit/mailman/message/27897393/
        Our thanks to Dr Michal Krompiec for the symmetrisation method and its implementation.
        Args:
            rdkit_mol:
                Molecule to find symmetry classes for (rdkit mol class object)
        return:
            A dict where the keys are the atom indices and the values are their types
            (type is arbitrarily based on index; only consistency is needed, no specific values)
        """

        # Check CIPRank is present for first atom (can assume it is present for all afterwards)
        if not rdkit_mol.GetAtomWithIdx(0).HasProp("_CIPRank"):
            Chem.AssignStereochemistry(
                rdkit_mol,
                cleanIt=True,
                force=True,
                flagPossibleStereoCenters=True,
            )

        # Array of ranks showing matching atoms
        cip_ranks = np.array(
            [int(atom.GetProp("_CIPRank")) for atom in rdkit_mol.GetAtoms()]
        )

        # Map the ranks to the atoms to produce a list of symmetrical atoms
        atom_symmetry_classes = [
            np.where(cip_ranks == rank)[0].tolist()
            for rank in range(max(cip_ranks) + 1)
        ]

        # Convert from list of classes to dict where each key is an atom and each value is its class (just a str)
        atom_symmetry_classes_dict = {}
        # i will be used to define the class (just index based)
        for i, sym_class in enumerate(atom_symmetry_classes):
            for atom in sym_class:
                atom_symmetry_classes_dict[atom] = str(i)

        return atom_symmetry_classes_dict

    @staticmethod
    def get_conformer_rmsd(
        rdkit_mol: Chem.Mol, ref_index: int, align_index: int
    ) -> float:
        """
        Get the rmsd between the current rdkit molecule and the coordinates provided.
        Args:
            rdkit_mol:
                rdkit representation of the molecule, conformer 0 is the base
            ref_index:
                The conformer index of the refernce
            align_index:
                the conformer index which should be aligned
        return:
            The RMSD value
        """

        return Chem.AllChem.GetConformerRMS(rdkit_mol, ref_index, align_index)

    @staticmethod
    def add_conformer(
        rdkit_mol: Chem.Mol, conformer_coordinates: np.ndarray
    ) -> Chem.Mol:
        """
        Add a new conformation to the rdkit molecule.
        Args:
            rdkit_mol:
                The rdkit molecule instance
            conformer_coordinates:
                A numpy array of the coordinates to be added
        return:
            The rdkit molecule with the conformer added
        """

        conformer = Chem.Conformer()
        for i, coord in enumerate(conformer_coordinates):
            atom_position = Point3D(*coord)
            conformer.SetAtomPosition(i, atom_position)

        rdkit_mol.AddConformer(conformer, assignId=True)

        return rdkit_mol


class ReadInput:
    """
    Called inside Ligand; used to handle reading any kind of input valid in QUBEKit
        QC JSON object
        SMILES string
        PDB, MOL2, XYZ file
    """

    def __init__(
        self,
        coords: Optional[np.ndarray] = None,
        rdkit_mol: Optional = None,
        name: Optional[str] = None,
    ):

        self.coords = coords
        self.rdkit_mol = rdkit_mol
        self.name = name

    @classmethod
    def from_smiles(
        cls, smiles: str, name: Optional[str] = None
    ) -> "ReadInput":
        """
        Make a ReadInput object which can be taken by the Ligand class to make the model.

        Note
        ----
        This method will generate a conformer for the molecule.

        Parameters
        ----------
        smiles:
            The smiles string which should be parsed by rdkit.
        name:
            The name that should be given to the molecule.
        """
        # Smiles string input
        rdkit_mol = RDKit.smiles_to_rdkit_mol(smiles_string=smiles, name=name)
        return cls(name=name, coords=None, rdkit_mol=rdkit_mol)

    @classmethod
    def from_file(cls, file_name: str) -> "ReadInput":
        """
        Read the input file using RDKit and return the molecule data.
        """
        input_file = Path(file_name)
        # if the file is not there raise an error
        if not input_file.exists():
            raise FileNotFoundError(
                f"{input_file.as_posix()} could not be found is this path correct?"
            )
        # xyz is a special case of file only internal readers catch
        if input_file.suffix == ".xyz":
            return cls.from_xyz(file_name=input_file.as_posix())
        # read the input with rdkit
        rdkit_mol = RDKit.file_to_rdkit_mol(file_path=input_file)
        return cls(
            rdkit_mol=rdkit_mol, coords=None, name=rdkit_mol.GetProp("_Name")
        )

    @classmethod
    def from_qc_json(cls, qc_json) -> "ReadInput":
        """
        Given a QC JSON object, extracts the topology, atoms and coords of the molecule.
        #TODO we need to be absle to read mapped smiles for this to work with stereochem and aromaticity
        """

        topology = nx.Graph()
        atoms = []

        for i, atom in enumerate(qc_json.symbols):
            atoms.append(
                Atom(
                    atomic_number=Element().number(atom),
                    atom_index=i,
                    atom_name=f"{atom}{i}",
                )
            )
            topology.add_node(i)

        for bond in qc_json.connectivity:
            topology.add_edge(*bond[:2])

        coords = (
            np.array(qc_json.geometry).reshape((len(atoms), 3)) * BOHR_TO_ANGS
        )
        return cls(name=None, rdkit_mol=None, coords=coords)

    @classmethod
    def from_xyz(cls, file_name: str) -> "ReadInput":
        """
        Internal xyz reader.
        Extracts the coords of the molecule.
        """

        traj_molecules = []
        coords = []

        with open(file_name) as xyz_file:
            lines = xyz_file.readlines()

            n_atoms = float(lines[0])

            for line in lines:
                line = line.split()
                # skip frame heading lines
                if len(line) <= 1 or "Iteration" in line:
                    continue

                coords.append([float(line[1]), float(line[2]), float(line[3])])

                if len(coords) == n_atoms:
                    # we have collected the molecule now store the frame
                    traj_molecules.append(np.array(coords))
                    coords = []

        coords = (
            traj_molecules[0] if len(traj_molecules) == 1 else traj_molecules
        )
        return cls(coords=coords, name=None, rdkit_mol=None)


class ReadInputProtein:
    """
    A class that specialises in reading Protein input files.
    #TODO are we better doing this with openmm or another tool?
    """

    def __init__(
        self,
        atoms: List[Atom],
        bonds: Optional[List[Bond]] = None,
        coords: Optional[np.ndarray] = None,
        residues: Optional[List[str]] = None,
        name: Optional[str] = None,
    ):
        self.atoms = atoms
        self.bonds = bonds
        self.coords = coords
        self.name = name
        self.residues = residues

    @classmethod
    def from_pdb(cls, file_name: str, name: Optional[str] = None):
        """
        Read the protein input pdb file.
        :return:
        """
        with open(file_name, "r") as pdb:
            lines = pdb.readlines()

        coords = []
        atoms = []
        bonds = []
        Residues = []

        # atom counter used for graph node generation
        atom_count = 0
        for line in lines:
            if "ATOM" in line or "HETATM" in line:
                atomic_symbol = str(line[76:78])
                atomic_symbol = re.sub("[0-9]+", "", atomic_symbol).strip()

                # If the element column is missing from the pdb, extract the atomic_symbol from the atom name.
                if not atomic_symbol:
                    atomic_symbol = str(line.split()[2])
                    atomic_symbol = re.sub("[0-9]+", "", atomic_symbol)

                # now make sure we have a valid element
                if (
                    atomic_symbol.lower() != "cl"
                    and atomic_symbol.lower() != "br"
                ):
                    atomic_symbol = atomic_symbol[0]

                atom_name = f"{atomic_symbol}{atom_count}"
                # TODO should we use a protein pdb package for this?
                qube_atom = Atom(
                    atomic_number=Element().number(atomic_symbol),
                    atom_index=atom_count,
                    atom_name=atom_name,
                    formal_charge=0,
                    aromatic=False,
                    bonds=[],
                )

                atoms.append(qube_atom)

                # also get the residue order from the pdb file so we can rewrite the file
                Residues.append(str(line.split()[3]))

                atom_count += 1
                coords.append(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ]
                )

            elif "CONECT" in line:
                conect_terms = line.split()
                for atom in conect_terms[2:]:
                    if int(atom):
                        bond = Bond(
                            atom1_index=int(conect_terms[1]) - 1,
                            atom2_index=int(atom) - 1,
                            bond_order=1,
                            aromatic=False,
                        )
                        bonds.append(bond)
                        atoms[int(conect_terms[1]) - 1].bonds.append(
                            int(atom) - 1
                        )

        coords = np.array(coords)
        residues = [res for res, group in groupby(Residues)]
        if name is None:
            name = Path(file_name).stem
        return cls(
            atoms=atoms,
            bonds=bonds,
            coords=coords,
            residues=residues,
            name=name,
        )


class Element:
    """
    Simple wrapper class for getting element info using RDKit.
    """

    @staticmethod
    def p_table() -> PeriodicTable:
        return GetPeriodicTable()

    @staticmethod
    def mass(identifier):
        pt = Element.p_table()
        return pt.GetAtomicWeight(identifier)

    @staticmethod
    def number(identifier):
        pt = Element.p_table()
        return pt.GetAtomicNumber(identifier)

    @staticmethod
    def name(identifier):
        pt = Element.p_table()
        return pt.GetElementSymbol(identifier)


class ReferenceData(BaseModel):
    """
    A basic reference data class.

    Note:
        Mixed units are used here due to the strange units in qdata.txt files.
    """

    geometry: Array[np.ndarray] = Field(
        ..., description="The geometry of the single point result in angstrom."
    )
    energy: float = Field(
        ...,
        description="The qm calculated energy for this geometry in hartree.",
    )
    gradient: Optional[Array[np.ndarray]] = Field(
        None,
        description="The gradient calculated at this geometry in hartree per bohr.",
    )

    class Config:
        validate_assignment = True
        json_encoders = {np.ndarray: lambda v: v.flatten().tolist()}

    @validator("geometry", "gradient")
    def _check_geom_grad(cls, array: np.ndarray) -> np.array:
        """
        Make sure that the geometry or gradient array is in the correct shape.
        """
        if array is None:
            return array
        else:
            return array.reshape((-1, 3))


class TorsionData(ReferenceData):
    """
    This is a reference data class for a grid point in a torsion drive scan.
    This class extends the normal data class with an angle attribute.
    """

    angle: int = Field(
        ..., description="The angle this data point was calculated at."
    )


class TorsionDriveData(BaseModel):
    """
    A container class to help store torsiondrive reference data the class is locked once made to ensure the data is valid.
    """

    grid_spacing: int = Field(
        15, description="The angle between grid points on a torsion drive."
    )
    torsion_drive_range: Tuple[int, int] = Field(
        (-165, 180),
        description="The torsion range this dihedral was scanned over.",
    )
    dihedral: Tuple[int, int, int, int] = Field(
        ...,
        description="The dihedral which was targeted during the torsion drive for fitting.",
    )
    reference_data: Dict[int, TorsionData] = Field(
        {},
        description="The list of single point torsion data points which detail the geometry, angle and energy.",
    )

    class Config:
        validate_assignment = True
        allow_mutation = False
        json_encoders = {np.ndarray: lambda v: v.flatten().tolist()}

    @validator("torsion_drive_range")
    def _check_scan_range(cls, scan_range: Tuple[int, int]) -> Tuple[int, int]:
        """
        Make sure the scan range is in order lowest to highest.
        """
        return tuple(sorted(scan_range))

    def to_file(self, file_name: str) -> None:
        """
        Write the object to file.
        """
        with open(file_name, "w") as output:
            output.write(self.json(indent=2))

    @property
    def max_angle(self) -> int:
        return self.torsion_drive_range[1]

    @property
    def min_angle(self) -> int:
        return self.torsion_drive_range[0]

    @property
    def possible_angles(self) -> List[int]:
        """
        Get all of the possible angle values for this grid spacing and scan range combination.
        """
        angles = [
            i
            for i in range(
                self.max_angle,
                self.min_angle - self.grid_spacing,
                -self.grid_spacing,
            )
            if i >= self.min_angle
        ]
        return angles

    @classmethod
    def from_qdata(
        cls,
        dihedral: Tuple[int, int, int, int],
        qdata_file: str = "qdata.txt",
        grid_spacing: int = 15,
        torsion_drive_range: Tuple[int, int] = (-165, 180),
    ) -> "TorsionDriveData":
        """
        Create a TorsionDrive Data class from a qdata.txt file and the target dihedral.

        Args:
            qdata_file: The path to the qdata file which should be read, be default this is qdata.txt.
            dihedral: The indices of the atoms which make up the target dihedral.
            grid_spacing: The spacing of the dihedral expected in the input file.
            torsion_drive_range: The torsion angle scan range expected in the input file.

        Returns:
            An instance of the class which holds the torsion drive data.
        """
        energies = []
        geometries = []
        angles = []
        # get the geometries and energies
        with open(qdata_file, "r") as qdata:
            for line in qdata.readlines():
                if "ENERGY" in line:
                    energies.append(float(line.split()[1]))
                elif "COORDS" in line:
                    coords = [float(x) for x in line.split()[1:]]
                    coords = np.array(coords).reshape((-1, 3))
                    # get the angle
                    angle = cls._measure_angle(
                        coords=coords, dihedral=dihedral
                    )
                    angles.append(angle)
                    geometries.append(coords)

        torsion_data = cls(
            grid_spacing=grid_spacing,
            torsion_drive_range=torsion_drive_range,
            dihedral=dihedral,
        )
        # now for each grid point at the reference data
        for i, angle in enumerate(angles):
            data_point = TorsionData(
                geometry=geometries[i], energy=energies[i], angle=angle
            )
            torsion_data.add_grid_point(grid_data=data_point)

        # make sure all of the reference data is consistent
        torsion_data.validate_angles()
        return torsion_data

    @staticmethod
    def _measure_angle(
        coords: np.ndarray, dihedral: Tuple[int, int, int, int]
    ) -> int:
        """
        For the given set of coords and dihedral atom indices measure the angle.

        Args:
            coords: The coordinates for which the dihedral should be measured in the shape (n, 3), where n is the number of atoms.
            dihedral: The indices of the atoms in the target dihedral.

        Note:
            The normal expected scan range for a torsiondrive is ~-170 to 180 so we flip a measured angle of -180 to 180.
        """
        # Calculate the dihedral angle in the molecule using the molecule data array.
        x1, x2, x3, x4 = [coords[atom_id] for atom_id in dihedral]
        b1, b2, b3 = x2 - x1, x3 - x2, x4 - x3
        t1 = np.linalg.norm(b2) * np.dot(b1, np.cross(b2, b3))
        t2 = np.dot(np.cross(b1, b2), np.cross(b2, b3))
        angle = round(np.degrees(np.arctan2(t1, t2)))
        # work around symmetry at 180
        if angle == -180:
            return 180
        else:
            return angle

    def validate_angles(self):
        """
        Make sure that the grid spacing and torsion scan range are consistent with the stored refernce data.
        """
        for i in self.possible_angles:
            if i not in self.reference_data:
                raise MissingReferenceData(
                    f"the torsion angle {i} has no reference data but is required for fitting."
                )

    def add_grid_point(self, grid_data: TorsionData) -> None:
        """
        Add some grid point data to the torsion drive dataset. Here we make sure an angle only appears once and is within
        the specified scan range.
        """
        if grid_data.angle in self.possible_angles:
            self.reference_data[grid_data.angle] = grid_data
        else:
            raise TorsionDriveDataError(
                f"Can not add data for torsion angle {grid_data.angle} as it is not consistent with a scan range off {self.torsion_drive_range} and grid spacing {self.grid_spacing}"
            )

    @property
    def central_bond(self) -> Tuple[int, int]:
        """
        Returns:
            The central bond tuple of the scanned torsion.
        """
        return tuple(self.dihedral[1:3])


class Molecule(SchemaBase):
    """Base class for ligands and proteins.

    The class is a simple representation of the molecule as a list of atom and bond objects, many attributes are then
    inferred from these core objects.
    """

    atoms: List[Atom] = Field(
        ...,
        description="A list of QUBEKit atom objects which make up the molecule.",
    )
    bonds: Optional[List[Bond]] = Field(
        None,
        description="The list of QUBEKit bond objects which connect the individual atoms.",
    )
    coordinates: Optional[Array[float]] = Field(
        None,
        description="A numpy arrary of the current cartesian positions of each atom in angstrom, this must be of size (n_atoms, 3)",
    )
    multiplicity: int = Field(
        1,
        description="The integer multiplicity of the molecule which is used in QM calculations.",
    )
    name: str = Field(
        "unk",
        description="An optional name string which will be used in all file IO calls by default.",
    )
    provenance: Dict[str, Any] = Field(
        dict(creator="QUBEKit", version=qubekit.__version__),
        description="Information on the version and method used to create this molecule.",
    )
    extra_sites: VirtualSiteGroup = Field(
        VirtualSiteGroup(),
        description="A force object which records any virtual sistes in the molecule",
    )
    BondForce: Union[HarmonicBondForce] = Field(
        HarmonicBondForce(),
        description="A force object which records bonded interactions between pairs of atoms",
    )
    AngleForce: Union[HarmonicAngleForce] = Field(
        HarmonicAngleForce(),
        description="A force object which records angle interactions between atom triplets.",
    )
    TorsionForce: PeriodicTorsionForce = Field(
        PeriodicTorsionForce(),
        description="A force object which records torsion interactions between atom quartets, using a periodic function.",
    )
    ImproperTorsionForce: PeriodicImproperTorsionForce = Field(
        PeriodicImproperTorsionForce(),
        description="A force group which records improper torsion interactions between atom quartets using a periodic function.",
    )
    RBTorsionForce: RBProperTorsionForce = Field(
        RBProperTorsionForce(),
        description="A force group which records torsion interactions between atom quartets using a RB function.",
    )
    ImproperRBTorsionForce: RBImproperTorsionForce = Field(
        RBImproperTorsionForce(),
        description="A force object which records improper torsion interactions between atom quartets using a RB function.",
    )
    NonbondedForce: Union[LennardJones126Force] = Field(
        LennardJones126Force(),
        description="A force group which records atom nonbonded parameters.",
    )
    chargemol_coords: Optional[Array[float]] = Field(
        None,
        description="The coordinates used to calculate the chargemol quantities, "
        "this is a reorientated conformation",
    )

    @validator("coordinates", "chargemol_coords")
    def _reshape_coords(
        cls, coordinates: Optional[np.ndarray]
    ) -> Optional[np.ndarray]:
        if coordinates is not None:
            return coordinates.reshape((-1, 3))
        else:
            return coordinates

    def __init__(
        self,
        atoms: List[Atom],
        bonds: Optional[List[Bond]] = None,
        coordinates: Optional[np.ndarray] = None,
        multiplicity: int = 1,
        name: str = "unk",
        routine: Optional[List[str]] = None,
        provenance: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        """
        Init the molecule using the basic information.

        Note:
            This method updates the provanance info.

        Args:
            atoms:
                A list of QUBEKit atom objects in the molecule.
            bonds:
                A list of QUBEKit bond objects in the molecule.
            coordinates:
                A numpy array of the current cartesian positions of each atom, this must be of size (n_atoms, 3)
            multiplicity:
                The integer multiplicity of the molecule which is used in QM calculations.
            name:
                An optional name string which will be used in all file IO calls by default.
            routine:
                The set of strings which encode the routine information used to create the molecule.

            bool; is the current execution starting from the beginning (False) or restarting (True)?
        """
        # the way the molecule was made
        method = routine or ["__init__"]
        if provenance is None:
            new_provenance = dict(
                creator="QUBEKit", version=qubekit.__version__, routine=method
            )
        else:
            # make sure we respect the provenance when parsing a json file
            new_provenance = provenance
            new_provenance["routine"] = new_provenance["routine"]

        super(Molecule, self).__init__(
            atoms=atoms,
            bonds=bonds,
            multiplicity=multiplicity,
            name=name,
            coordinates=coordinates,
            provenance=new_provenance,
            **kwargs,
        )

    def to_topology(self) -> nx.Graph:
        """
        Build a networkx representation of the molecule.
        TODO add other attributes to the graph?
        """
        graph = nx.Graph()
        for atom in self.atoms:
            graph.add_node(atom.atom_index)

        for bond in self.bonds:
            graph.add_edge(bond.atom1_index, bond.atom2_index)
        return graph

    def to_file(self, file_name: str) -> None:
        """
        Write the molecule object to file working out the file type from the extension.
        Works with PDB, MOL, SDF, XYZ, Json any other we want?
        """
        if ".json" in file_name:
            with open(file_name, "w") as output:
                output.write(self.json(indent=2))
        else:
            return RDKit.mol_to_file(
                rdkit_mol=self.to_rdkit(), file_name=file_name
            )

    def generate_conformers(self, n_conformers: int) -> List[np.ndarray]:
        """
        Generate a list of conformers and return them this includes the input conformer

        Args:
            n_conformers: The number of conformers which should be generated.
        """
        rd_mol = self.to_rdkit()
        return RDKit.generate_conformers(
            rdkit_mol=rd_mol, conformer_no=n_conformers
        )

    def to_multiconformer_file(
        self, file_name: str, positions: List[np.ndarray]
    ) -> None:
        """
        Write the molecule to a file allowing multipule conformers.

        As the ligand object only holds one set of coordinates at once a list of coords can be passed here to allow
        multiconformer support.

        Args:
            file_name:
                The name of the file that should be created, the type is inferred by the suffix.
            positions:
                A list of Cartesian coordinates of shape (n_atoms, 3).
        """
        rd_mol = self.to_rdkit()
        rd_mol.RemoveAllConformers()
        # add the conformers
        if not isinstance(positions, list):
            positions = [
                positions,
            ]
        for conformer in positions:
            RDKit.add_conformer(
                rdkit_mol=rd_mol, conformer_coordinates=conformer
            )
        return RDKit.mol_to_multiconformer_file(
            rdkit_mol=rd_mol, file_name=file_name
        )

    def get_atom_with_name(self, name):
        """
        Search through the molecule for an atom with that name and return it when found
        :param name: The name of the atom we are looking for
        :return: The QUBE Atom object with the name
        """

        for atom in self.atoms:
            if atom.atom_name == name:
                return atom
        raise AttributeError("No atom found with that name.")

    def get_bond_between(self, atom1_index: int, atom2_index: int) -> Bond:
        """
        Try and find a bond between the two atom indices.

        The bond may not have the atoms in the expected order.

        Args:
            atom1_index:
                The index of the first atom in the atoms list.
            atom2_index:
                The index of the second atom in the atoms list.

        Returns:
            The bond object between the two target atoms.

        Raises:
            TopologyMismatch:
                When no bond can be found between the atoms.
        """
        target = [atom1_index, atom2_index]
        for bond in self.bonds:
            if bond.atom1_index in target and bond.atom2_index in target:
                return bond
        raise TopologyMismatch(
            f"There is no bond between atoms {atom1_index} and {atom2_index} in this molecule."
        )

    @property
    def has_unique_atom_names(self) -> bool:
        """
        Check if the molecule has unique atom names or not this will help with pdb file writing.
        """
        atom_names = set([atom.atom_name for atom in self.atoms])
        if len(atom_names) == self.n_atoms:
            return True
        return False

    @property
    def improper_torsions(self) -> Optional[List[Tuple[int, int, int, int]]]:
        """A list of improper atom tuples where the first atom is central."""

        improper_torsions = []
        topology = self.to_topology()
        for node in topology.nodes:
            near = sorted(list(nx.neighbors(topology, node)))
            # if the atom has 3 bonds it could be an improper
            # Check if an sp2 carbon or N
            if len(near) == 3 and (
                self.atoms[node].atomic_symbol == "C"
                or self.atoms[node].atomic_symbol == "N"
            ):
                # Store each combination of the improper torsion
                improper_torsions.append((node, near[0], near[1], near[2]))
        return improper_torsions or None

    @property
    def n_improper_torsions(self) -> int:
        """The number of unique improper torsions."""
        impropers = self.improper_torsions
        if impropers is None:
            return 0
        return len(impropers)

    @property
    def angles(self) -> Optional[List[Tuple[int, int, int]]]:
        """A List of angles from the topology."""

        angles = []
        topology = self.to_topology()
        for node in topology.nodes:
            bonded = sorted(list(nx.neighbors(topology, node)))

            # Check that the atom has more than one bond
            if len(bonded) < 2:
                continue

            # Find all possible angle combinations from the list
            for i in range(len(bonded)):
                for j in range(i + 1, len(bonded)):
                    atom1, atom3 = bonded[i], bonded[j]
                    angles.append((atom1, node, atom3))
        return angles or None

    @property
    def charge(self) -> int:
        """
        Return the integer charge of the molecule as the sum of the formal charge.
        """
        return sum([atom.formal_charge for atom in self.atoms])

    @property
    def n_angles(self) -> int:
        """The number of angles in the molecule."""
        angles = self.angles
        if angles is None:
            return 0
        return len(angles)

    def measure_bonds(self) -> Dict[Tuple[int, int], float]:
        """
        Find the length of all bonds in the molecule for the given conformer in  angstroms.

        Returns:
            A dictionary of the bond lengths stored by bond tuple.
        """

        bond_lengths = {}

        for bond in self.bonds:
            atom1 = self.coordinates[bond.atom1_index]
            atom2 = self.coordinates[bond.atom2_index]
            edge = (bond.atom1_index, bond.atom2_index)
            bond_lengths[edge] = np.linalg.norm(atom2 - atom1)

        return bond_lengths

    @property
    def n_bonds(self) -> int:
        """The number of bonds in the topology."""
        bonds = self.bonds
        if bonds is None:
            return 0
        return len(bonds)

    @property
    def dihedrals(
        self,
    ) -> Optional[Dict[Tuple[int, int], List[Tuple[int, int, int, int]]]]:
        """A list of all possible dihedrals that can be found in the topology."""

        dihedrals = {}
        topology = self.to_topology()
        # Work through the network using each edge as a central dihedral bond
        for edge in topology.edges:

            for start in list(nx.neighbors(topology, edge[0])):

                # Check atom not in main bond
                if start != edge[0] and start != edge[1]:

                    for end in list(nx.neighbors(topology, edge[1])):

                        # Check atom not in main bond
                        if end != edge[0] and end != edge[1]:

                            if edge not in dihedrals:
                                # Add the central edge as a key the first time it is used
                                dihedrals[edge] = [
                                    (start, edge[0], edge[1], end)
                                ]

                            else:
                                # Add the tuple to the correct key.
                                dihedrals[edge].append(
                                    (start, edge[0], edge[1], end)
                                )

        return dihedrals or None

    @property
    def n_dihedrals(self) -> int:
        """The total number of dihedrals in the molecule."""
        dihedrals = self.dihedrals
        if dihedrals is None:
            return 0
        return sum([len(torsions) for torsions in dihedrals.values()])

    def find_rotatable_bonds(
        self, smirks_to_remove: Optional[List[str]] = None
    ) -> Optional[List[Bond]]:
        """
        Args:
            smirks_to_remove:
                Optional list of smirks patterns which will be discarded
                from the rotatable bonds
        Find all rotatable bonds in the molecule.
        Remove any groups which are not relevant for torsion scans.
            e.g. methyl / amine groups
        return:
            The rotatable bonds in the molecule to be used for torsion scans.
        """

        rotatable_bond_smarts = "[!$(*#*)&!D1:1]-&!@[!$(*#*)&!D1:2]"

        rotatable_matches = self.get_smarts_matches(rotatable_bond_smarts)
        if rotatable_matches is None:
            return None

        if smirks_to_remove is not None:
            for smirk in smirks_to_remove:
                matches_to_remove = self.get_smarts_matches(smirk)
                if matches_to_remove is not None:
                    for match in matches_to_remove:
                        try:
                            rotatable_matches.remove(match)
                        except ValueError:
                            try:
                                # If the match is not in the list, it may be in backwards
                                rotatable_matches.remove(
                                    tuple(reversed(match))
                                )
                            except ValueError:
                                continue

        # gather a list of bond instances to return
        rotatable_bonds = [
            self.get_bond_between(*bond) for bond in rotatable_matches
        ]

        return rotatable_bonds or None

    @property
    def n_rotatable_bonds(self) -> int:
        """The number of rotatable bonds."""
        rotatable_bonds = self.find_rotatable_bonds()
        if rotatable_bonds is None:
            return 0
        return len(rotatable_bonds)

    def symmetrise_nonbonded_parameters(self) -> bool:
        """
        Symmetrise all non-bonded force group parameters.

        Using the CIP rankings from RDKit apply symmetry to the non-bonded force group.

        Important:
            We respect the predefined parameters in the non-bonded force group which can be symmetrised.
        """
        # group atom types as they are in a different format to other types
        atom_types = {}
        for atom_index, cip_type in self.atom_types.items():
            atom_types.setdefault(cip_type, []).append((atom_index,))
        for atoms in atom_types.items():
            self._symmetrise_parameters(
                force_group=self.NonbondedForce, parameter_keys=atoms
            )

        return True

    def symmetrise_bonded_parameters(self) -> bool:
        """
        Symmetrise all bond and angle force group parameters.

        Using the CIP rankings from RDKit apply symmetry to the bond and angle force groups.

        Important:
            We respect the predefined parameters in the bond/angle force group which can be symmetrised.
        """

        for bonds in self.bond_types.values():
            self._symmetrise_parameters(
                force_group=self.BondForce, parameter_keys=bonds
            )

        for angles in self.angle_types.values():
            self._symmetrise_parameters(
                force_group=self.AngleForce, parameter_keys=angles
            )

        return True

    def _symmetrise_parameters(
        self,
        force_group: BaseForceGroup,
        parameter_keys: List[Tuple[int, ...]],
    ):
        """
        Internal method which applies symmetry to a group of parameter references in a particular force group.

        Args:
            force_group: The force group we should query for parameters.
            parameter_keys: The list of atom indices tuples that the symmetry should be applied to.
        """

        symmetry_attrs = force_group.symmetry_parameters()

        raw_parameter_values = {}
        for parameter_key in parameter_keys:
            param = force_group[parameter_key]
            for attr in symmetry_attrs:
                raw_parameter_values.setdefault(attr, []).append(
                    getattr(param, attr)
                )

        # now average the raw values
        for key, value in raw_parameter_values.items():
            raw_parameter_values[key] = np.array(value).mean()

        # now set back
        for parameter_key in parameter_keys:
            force_group.create_parameter(
                atoms=parameter_key, **raw_parameter_values
            )

    def measure_dihedrals(
        self,
    ) -> Optional[Dict[Tuple[int, int, int, int], float]]:
        """
        For the given conformation measure the dihedrals in the topology in degrees.
        """
        dihedrals = self.dihedrals
        if dihedrals is None:
            return None

        dih_phis = {}

        for val in dihedrals.values():
            for torsion in val:
                # Calculate the dihedral angle in the molecule using the molecule data array.
                x1, x2, x3, x4 = [
                    self.coordinates[torsion[i]] for i in range(4)
                ]
                b1, b2, b3 = x2 - x1, x3 - x2, x4 - x3
                t1 = np.linalg.norm(b2) * np.dot(b1, np.cross(b2, b3))
                t2 = np.dot(np.cross(b1, b2), np.cross(b2, b3))
                dih_phis[torsion] = np.degrees(np.arctan2(t1, t2))

        return dih_phis

    def measure_angles(self) -> Optional[Dict[Tuple[int, int, int], float]]:
        """
        For the given conformation measure the angles in the topology in degrees.
        """
        angles = self.angles
        if angles is None:
            return None

        angle_values = {}

        for angle in angles:
            x1 = self.coordinates[angle[0]]
            x2 = self.coordinates[angle[1]]
            x3 = self.coordinates[angle[2]]
            b1, b2 = x1 - x2, x3 - x2
            cosine_angle = np.dot(b1, b2) / (
                np.linalg.norm(b1) * np.linalg.norm(b2)
            )
            angle_values[angle] = np.degrees(np.arccos(cosine_angle))

        return angle_values

    @property
    def n_atoms(self) -> int:
        """
        Calculate the number of atoms.
        """
        return len(self.atoms)

    def write_parameters(self, file_name: str):
        """
        Take the molecule's parameter set and write an xml file for the molecule.
        """

        tree = self._build_forcefield().getroot()
        messy = ET.tostring(tree, "utf-8")

        pretty_xml_as_string = parseString(messy).toprettyxml(indent="")

        with open(file_name, "w") as xml_doc:
            xml_doc.write(pretty_xml_as_string)

    def _build_forcefield(self):
        """
        Separates the parameters and builds an xml tree ready to be used.

        TODO how do we support OPLS combination rules.
        Important:
            The ordering here should not be changed due to the way sites have to be added.
        """

        # Create XML layout
        root = ET.Element("ForceField")

        ET.SubElement(
            root,
            "QUBEKit",
            attrib={
                "Version": qubekit.__version__,
                "Date": datetime.now().strftime("%Y_%m_%d"),
            },
        )
        AtomTypes = ET.SubElement(root, "AtomTypes")
        Residues = ET.SubElement(root, "Residues")

        resname = "QUP" if self.__class__.__name__ == "Protein" else "MOL"
        Residue = ET.SubElement(Residues, "Residue", name=resname)
        # declare atom `types` and properties
        for atom in self.atoms:
            atom_type = f"QUBE_{atom.atom_index}"
            ET.SubElement(
                AtomTypes,
                "Type",
                attrib={
                    "name": atom_type,
                    "class": str(atom.atom_index),
                    "element": atom.atomic_symbol,
                    "mass": str(atom.atomic_mass),
                },
            )

            ET.SubElement(
                Residue,
                "Atom",
                attrib={"name": atom.atom_name, "type": atom_type},
            )

        # add sites to Atomtypes, topology and nonbonded
        for i, site in enumerate(self.extra_sites, start=1):
            site_name = f"v-site{i}"
            site_class = f"X{i}"
            ET.SubElement(
                AtomTypes,
                "Type",
                attrib={"name": site_name, "class": site_class, "mass": "0"},
            )
            # for some reason we swap name and class here but it works !
            ET.SubElement(
                Residue, "Atom", attrib={"name": site_class, "type": site_name}
            )

        BondForce = ET.SubElement(
            root,
            self.BondForce.openmm_group(),
            attrib=self.BondForce.xml_data(),
        )
        for parameter in self.BondForce:
            ET.SubElement(
                BondForce, parameter.openmm_type(), attrib=parameter.xml_data()
            )
            ET.SubElement(
                Residue,
                "Bond",
                attrib={
                    "from": str(parameter.atoms[0]),
                    "to": str(parameter.atoms[1]),
                },
            )
        AngleForce = ET.SubElement(
            root,
            self.AngleForce.openmm_group(),
            attrib=self.AngleForce.xml_data(),
        )
        for parameter in self.AngleForce:
            ET.SubElement(
                AngleForce,
                parameter.openmm_type(),
                attrib=parameter.xml_data(),
            )
        if (
            self.TorsionForce.n_parameters > 0
            or self.ImproperTorsionForce.n_parameters > 0
        ):
            TorsionForce = ET.SubElement(
                root,
                self.TorsionForce.openmm_group(),
                attrib=self.TorsionForce.xml_data(),
            )
            for parameter in self.TorsionForce:
                ET.SubElement(
                    TorsionForce,
                    parameter.openmm_type(),
                    attrib=parameter.xml_data(),
                )
            for parameter in self.ImproperTorsionForce:
                ET.SubElement(
                    TorsionForce,
                    parameter.openmm_type(),
                    attrib=parameter.xml_data(),
                )
        if (
            self.RBTorsionForce.n_parameters > 0
            or self.ImproperRBTorsionForce.n_parameters > 0
        ):
            RBTorsion = ET.SubElement(
                root,
                self.RBTorsionForce.openmm_group(),
                attrib=self.RBTorsionForce.xml_data(),
            )
            for parameter in self.RBTorsionForce:
                ET.SubElement(
                    RBTorsion,
                    parameter.openmm_type(),
                    attrib=parameter.xml_data(),
                )
            for parameter in self.ImproperRBTorsionForce:
                ET.SubElement(
                    RBTorsion,
                    parameter.openmm_type(),
                    attrib=parameter.xml_data(),
                )

        # now we add more site info after general bonding
        for i, site in enumerate(self.extra_sites):
            site_data = site.xml_data()
            # we have to add its global index
            site_data["index"] = str(i + self.n_atoms)
            ET.SubElement(Residue, site.openmm_type(), attrib=site_data)

        NonbondedForce = ET.SubElement(
            root,
            self.NonbondedForce.openmm_group(),
            attrib=self.NonbondedForce.xml_data(),
        )
        for parameter in self.NonbondedForce:
            ET.SubElement(
                NonbondedForce,
                parameter.openmm_type(),
                attrib=parameter.xml_data(),
            )

        for i, site in enumerate(self.extra_sites, start=1):
            site_name = f"v-site{i}"
            ET.SubElement(
                NonbondedForce,
                "Atom",
                attrib={
                    "charge": str(site.charge),
                    "epsilon": "0",
                    "sigma": "1",
                    "type": site_name,
                },
            )

        return ET.ElementTree(root)

    @property
    def bond_types(self) -> Dict[str, List[Tuple[int, int]]]:
        """
        Using the symmetry dict, give each bond a code. If any codes match, the bonds can be symmetrised.
        e.g. bond_symmetry_classes = {(0, 3): '2-0', (0, 4): '2-0', (0, 5): '2-0' ...}
        all of the above bonds (tuples) are of the same type (methyl H-C bonds in same region)
        This dict is then used to produce bond_types.
        bond_types is just a dict where the keys are the string code from above and the values are all
        of the bonds with that particular type.
        """
        atom_types = self.atom_types
        bond_symmetry_classes = {}
        for bond in self.bonds:
            bond_symmetry_classes[(bond.atom1_index, bond.atom2_index)] = (
                f"{atom_types[bond.atom1_index]}-"
                f"{atom_types[bond.atom2_index]}"
            )

        bond_types = {}
        for key, val in bond_symmetry_classes.items():
            bond_types.setdefault(val, []).append(key)

        bond_types = self._cluster_types(bond_types)
        return bond_types

    @property
    def angle_types(self) -> Dict[str, List[Tuple[int, int, int]]]:
        """
        Using the symmetry dict, give each angle a code. If any codes match, the angles can be symmetrised.
        e.g. angle_symmetry_classes = {(1, 0, 3): '3-2-0', (1, 0, 4): '3-2-0', (1, 0, 5): '3-2-0' ...}
        all of the above angles (tuples) are of the same type (methyl H-C-H angles in same region)
        angle_types is just a dict where the keys are the string code from the above and the values are all
        of the angles with that particular type.
        """
        atom_types = self.atom_types
        angle_symmetry_classes = {}
        for angle in self.angles:
            angle_symmetry_classes[angle] = (
                f"{atom_types[angle[0]]}-"
                f"{atom_types[angle[1]]}-"
                f"{atom_types[angle[2]]}"
            )

        angle_types = {}
        for key, val in angle_symmetry_classes.items():
            angle_types.setdefault(val, []).append(key)

        angle_types = self._cluster_types(angle_types)
        return angle_types

    @property
    def dihedral_types(self) -> Dict[str, List[Tuple[int, int, int, int]]]:
        """
        Using the symmetry dict, give each dihedral a code. If any codes match, the dihedrals can be clustered and their
        parameters should be the same, this is to be used in dihedral fitting so all symmetry equivalent dihedrals are
        optimised at the same time. dihedral_equiv_classes = {(0, 1, 2 ,3): '1-1-2-1'...} all of the tuples are the
        dihedrals index by topology and the strings are the symmetry equivalent atom combinations.
        """
        atom_types = self.atom_types
        dihedral_symmetry_classes = {}
        for dihedral_set in self.dihedrals.values():
            for dihedral in dihedral_set:
                dihedral_symmetry_classes[tuple(dihedral)] = (
                    f"{atom_types[dihedral[0]]}-"
                    f"{atom_types[dihedral[1]]}-"
                    f"{atom_types[dihedral[2]]}-"
                    f"{atom_types[dihedral[3]]}"
                )

        dihedral_types = {}
        for key, val in dihedral_symmetry_classes.items():
            dihedral_types.setdefault(val, []).append(key)

        dihedral_types = self._cluster_types(dihedral_types)
        return dihedral_types

    @property
    def improper_types(self) -> Dict[str, List[Tuple[int, int, int, int]]]:
        """Using the atom symmetry types work out the improper types."""

        atom_types = self.atom_types
        improper_symmetry_classes = {}
        for dihedral in self.improper_torsions:
            improper_symmetry_classes[tuple(dihedral)] = (
                f"{atom_types[dihedral[0]]}-"
                f"{atom_types[dihedral[1]]}-"
                f"{atom_types[dihedral[2]]}-"
                f"{atom_types[dihedral[3]]}"
            )

        improper_types = {}
        for key, val in improper_symmetry_classes.items():
            improper_types.setdefault(val, []).append(key)

        improper_types = self._cluster_types(improper_types)
        return improper_types

    @staticmethod
    def _cluster_types(equiv_classes):
        """
        Function that helps the bond angle and dihedral class finders in clustering the types based on the forward and
        backward type strings.
        :return: clustered equiv class
        """

        new_classes = {}
        for key, item in equiv_classes.items():
            try:
                new_classes[key].extend(item)
            except KeyError:
                try:
                    new_classes[key[::-1]].extend(item)
                except KeyError:
                    new_classes[key] = item

        return new_classes

    @property
    def atom_types(self) -> Dict[int, str]:
        """Returns a dictionary of atom indices mapped to their class or None if there is no rdkit molecule."""

        return RDKit.find_symmetry_classes(self.to_rdkit())

    def to_rdkit(self) -> Chem.Mol:
        """
        Generate an rdkit representation of the QUBEKit ligand object.

        Here we build the molecule and assign the stereochemistry using the coordinates as we should always have a set of coordinates in the model.
        This allows us to skip complicated local vs global stereo chemistry checks however this could break in future.

        Returns:
            An rdkit representation of the molecule.
        """
        from qubekit.utils.helpers import _assert_wrapper

        # TODO what properties should be put in the rdkit molecule? Multiplicity?
        # make an editable molecule
        rd_mol = Chem.RWMol()
        if self.name is not None:
            rd_mol.SetProp("_Name", self.name)

        # when building the molecule we have to loop multiple times
        # so always make sure the indexing is the same in qube and rdkit
        for atom in self.atoms:
            rd_index = rd_mol.AddAtom(atom.to_rdkit())
            assert rd_index == atom.atom_index

        # now we need to add each bond, can not make a bond from python currently
        for bond in self.bonds:
            rd_mol.AddBond(*bond.indices)
            # now get the bond back to edit it
            rd_bond: Chem.Bond = rd_mol.GetBondBetweenAtoms(*bond.indices)
            rd_bond.SetIsAromatic(bond.aromatic)
            rd_bond.SetBondType(bond.rdkit_type)

        Chem.SanitizeMol(
            rd_mol,
            Chem.SANITIZE_ALL
            ^ Chem.SANITIZE_ADJUSTHS
            ^ Chem.SANITIZE_SETAROMATICITY,
        )
        # must use openff MDL model for compatibility
        Chem.SetAromaticity(rd_mol, Chem.AromaticityModel.AROMATICITY_MDL)

        # conformers
        rd_mol = RDKit.add_conformer(
            rdkit_mol=rd_mol, conformer_coordinates=self.coordinates
        )
        Chem.AssignStereochemistryFrom3D(rd_mol)

        # now we should check that the stereo has not been broken
        for rd_atom in rd_mol.GetAtoms():
            index = rd_atom.GetIdx()
            qb_atom = self.atoms[index]
            if qb_atom.stereochemistry is not None:
                with _assert_wrapper(StereoChemistryError):
                    assert qb_atom.stereochemistry == rd_atom.GetProp(
                        "_CIPCode"
                    ), f"StereoChemistry incorrect expected {qb_atom.stereochemistry} got {rd_atom.GetProp('_CIPCode')} for atom {qb_atom}"

        for rd_bond in rd_mol.GetBonds():
            index = rd_bond.GetIdx()
            qb_bond = self.bonds[index]
            if qb_bond.stereochemistry is not None:
                rd_bond.SetStereo(qb_bond.rdkit_stereo)
            rd_stereo = rd_bond.GetStereo()
            if qb_bond.stereochemistry == "E":
                with _assert_wrapper(StereoChemistryError):
                    assert (
                        rd_stereo == Chem.BondStereo.STEREOE
                    ), f"StereoChemistry incorrect expected E got {rd_stereo}"
            elif qb_bond.stereochemistry == "Z":
                with _assert_wrapper(StereoChemistryError):
                    assert (
                        rd_stereo == Chem.BondStereo.STEREOZ
                    ), f"StereoChemistry incorrect expected Z got {rd_stereo}"

        return Chem.Mol(rd_mol)

    def get_smarts_matches(
        self, smirks: str
    ) -> Optional[List[Tuple[int, ...]]]:
        """
        Get substructure matches for a mapped SMARTS pattern.

        Args:
            smirks:
                The mapped SMARTS pattern that should be used to query the molecule.

        Returns:
            `None` if there are no matches, else a list of tuples of atom indices which match the tagged atoms in
            the SMARTS pattern. These are returned in the same order.
        """
        matches = RDKit.get_smirks_matches(
            rdkit_mol=self.to_rdkit(), smirks=smirks
        )
        if not matches:
            return None
        return matches

    def add_qm_scan(self, scan_data: TorsionDriveData) -> None:
        """
        Save the torsion drive data into the ligand object.
        """
        if scan_data.__class__ != TorsionDriveData:
            raise MissingReferenceData(
                "The reference data must be in the form of the torsion drive data class."
            )
        else:
            if self.qm_scans is None:
                self.qm_scans = []
            self.qm_scans.append(scan_data)

    def openmm_coordinates(self) -> unit.Quantity:
        """
        Convert the coordinates to an openMM quantity.

        Build a single set of coordinates for the molecule that work in openMM.
        Note this must be a single conformer, if multiple are given only the first is used.

        Returns:
            A openMM quantity wrapped array of the coordinates in angstrom.
        """
        return unit.Quantity(self.coordinates, unit.angstroms)

    def fix_net_charge(self):
        """
        Ensure the total is exactly equal to the ideal net charge of the molecule.
        If net charge is not an integer value, MM simulations can (ex/im)plode.
        """

        decimal.setcontext(decimal.Context(prec=7))
        round_to = decimal.Decimal(10) ** -6

        for param in self.NonbondedForce:
            param.charge = param.charge.quantize(round_to)

        atom_charges = sum(param.charge for param in self.NonbondedForce)
        extra = self.charge - atom_charges

        if self.extra_sites is not None:
            for site in self.extra_sites:
                site.charge = site.charge.quantize(round_to)
                extra -= site.charge

        if extra:
            last_atom_index = self.n_atoms - 1
            self.NonbondedForce[(last_atom_index,)].charge += extra


class Ligand(Molecule):
    """
    The Ligand class seperats from protiens as we add fields to store QM calculations, such as the hessian and add more
    rdkit support methods.
    """

    hessian: Optional[Array[float]] = Field(
        None,
        description="The hessian matrix calculated for this molecule at the QM optimised geometry.",
    )
    qm_scans: Optional[List[TorsionDriveData]] = Field(
        None,
        description="The list of reference torsiondrive results which we can fit against.",
    )
    wbo: Optional[Array[float]] = Field(
        None,
        description="The WBO matrix calculated at the QM optimised geometry.",
    )

    def __init__(
        self,
        atoms: List[Atom],
        bonds: Optional[List[Bond]] = None,
        coordinates: Optional[np.ndarray] = None,
        multiplicity: int = 1,
        name: str = "unk",
        routine: Optional[Set] = None,
        provenance: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        super(Ligand, self).__init__(
            atoms=atoms,
            bonds=bonds,
            coordinates=coordinates,
            multiplicity=multiplicity,
            name=name,
            routine=routine,
            provenance=provenance,
            **kwargs,
        )
        # make sure we have unique atom names
        self._validate_atom_names()

    @validator("hessian", "wbo", allow_reuse=True)
    def _reshape_matrix(
        cls, matrix: Optional[np.ndarray]
    ) -> Optional[np.ndarray]:
        if matrix is not None:
            if len(matrix.shape) == 1:
                # the matrix is a flat list
                # so we need to make the matrix to be square
                length = int(np.sqrt(matrix.shape[0]))
                return matrix.reshape((length, length))
        return matrix

    @classmethod
    def from_rdkit(
        cls,
        rdkit_mol: Chem.Mol,
        name: Optional[str] = None,
        multiplicity: int = 1,
    ) -> "Ligand":
        """
        Build an instance of a qubekit ligand directly from an rdkit molecule.

        Args:
            rdkit_mol:
                An instance of an rdkit.Chem.Mol from which the QUBEKit ligand should be built.
            name:
                The name that should be assigned to the molecule, this will overwrite any name already assigned.
            multiplicity:
                The multiplicity of the molecule, used in QM calculations.
        """
        if name is None:
            if rdkit_mol.HasProp("_Name"):
                name = rdkit_mol.GetProp("_Name")

        atoms = []
        bonds = []
        # Collect the atom names and bonds
        for rd_atom in rdkit_mol.GetAtoms():
            # make and atom
            qb_atom = Atom.from_rdkit(rd_atom=rd_atom)
            atoms.append(qb_atom)

        # now we need to make a list of bonds
        for rd_bond in rdkit_mol.GetBonds():
            qb_bond = Bond.from_rdkit(rd_bond=rd_bond)
            bonds.append(qb_bond)

        coords = rdkit_mol.GetConformer().GetPositions()
        bonds = bonds or None
        # method use to make the molecule
        routine = ["QUBEKit.ligand.from_rdkit"]
        return cls(
            atoms=atoms,
            bonds=bonds,
            coordinates=coords,
            multiplicity=multiplicity,
            name=name,
            routine=routine,
        )

    def has_ub_terms(self) -> bool:
        """Return `True` if the molecule has Ure-Bradly terms, as there are forces between non-bonded atoms."""
        for bond in self.BondForce:
            try:
                self.get_bond_between(*bond.atoms)
            except TopologyMismatch:
                return True
        return False

    def _to_ub_pdb(self, file_name: Optional[str] = None) -> None:
        """A privet method to write the molecule to a non-standard pdb file with connections for Urey-Bradly terms."""
        openmm_top = self.to_openmm_topology()
        PDBFile.writeFile(
            topology=openmm_top,
            positions=self.openmm_coordinates(),
            file=open(f"{file_name or self.name}.pdb", "w"),
        )

    @staticmethod
    def _check_file_name(file_name: str) -> None:
        """
        Make sure that if an unsupported file type is passed we can not make a molecule from it.
        """
        if ".xyz" in file_name:
            raise FileTypeError(
                "XYZ files can not be used to build ligands due to ambiguous bonding, "
                "please use pdb, mol, mol2 or smiles as input."
            )

    @classmethod
    def from_file(cls, file_name: str, multiplicity: int = 1) -> "Ligand":
        """
        Build a ligand from a supported input file.

        Args:
            file_name:
                The abs path to the file including the extension which determines how the file is read.
            multiplicity:
                The multiplicity of the molecule which is required for QM calculations.
        """
        cls._check_file_name(file_name=file_name)
        input_data = ReadInput.from_file(file_name=file_name)
        ligand = cls.from_rdkit(
            rdkit_mol=input_data.rdkit_mol,
            name=input_data.name,
            multiplicity=multiplicity,
        )
        # now edit the routine to include this call
        ligand.provenance["routine"].extend(
            ["QUBEKit.ligand.from_file", os.path.abspath(file_name)]
        )
        return ligand

    @classmethod
    def from_smiles(
        cls, smiles_string: str, name: str, multiplicity: int = 1
    ) -> "Ligand":
        """
        Build the ligand molecule directly from a non mapped smiles string.

        Args:
            smiles_string:
                The smiles string from which a molecule instance should be made.
            name:
                The name that should be assigned to the molecule.
            multiplicity:
                The multiplicity of the molecule, important for QM calculations.
        """
        input_data = ReadInput.from_smiles(smiles=smiles_string, name=name)
        ligand = cls.from_rdkit(
            rdkit_mol=input_data.rdkit_mol,
            name=name,
            multiplicity=multiplicity,
        )
        # now edit the routine to include this command
        ligand.provenance["routine"].extend(
            ["QUBEKit.ligand.from_smiles", smiles_string]
        )
        return ligand

    def to_openmm_topology(self) -> Topology:
        """
        Convert the Molecule to a OpenMM topology representation.

        We assume we have a single molecule so a single chain is made with a single residue.
        Note this will not work with proteins as we will need to have distinct residues.

        Returns:
            An openMM topology object which can be used to construct a system.
        """
        topology = Topology()
        bond_types = {1: Single, 2: Double, 3: Triple}
        chain = topology.addChain()
        # create a molecule specific residue
        residue = topology.addResidue(name=self.name, chain=chain)
        # add atoms and keep track so we can add bonds
        top_atoms = []
        for atom in self.atoms:
            element = Element.getByAtomicNumber(atom.atomic_number)
            top_atom = topology.addAtom(
                name=atom.atom_name, element=element, residue=residue
            )
            top_atoms.append(top_atom)
        for bond in self.bonds:
            atom1 = top_atoms[bond.atom1_index]
            atom2 = top_atoms[bond.atom2_index]
            # work out the type
            if bond.aromatic:
                b_type = Aromatic
            else:
                b_type = bond_types[bond.bond_order]
            topology.addBond(
                atom1=atom1, atom2=atom2, type=b_type, order=bond.bond_order
            )
        # now check for Urey-Bradley terms
        for bond in self.BondForce:
            try:
                self.get_bond_between(*bond.atoms)
            except TopologyMismatch:
                # the bond is not in the bond list so add it as a u-b term
                atom1 = top_atoms[bond.atoms[0]]
                atom2 = top_atoms[bond.atoms[1]]
                # this is a fake bond used for U-B terms.
                topology.addBond(
                    atom1=atom1, atom2=atom2, type=Single, order=1
                )

        return topology

    def to_smiles(
        self,
        isomeric: bool = True,
        explicit_hydrogens: bool = True,
        mapped: bool = False,
    ) -> str:
        """
        Create a canonical smiles representation for the molecule based on the input settings.

        Args:
            isomeric:
                If the smiles string should encode stereochemistry `True` or not `False`.
            explicit_hydrogens:
                If hydrogens should be explicitly encoded into the smiles string `True` or not `False`.
            mapped:
                If the smiles should encode the original atom ordering `True` or not `False` as this might be different
                from the canonical ordering.

        Returns:
            A smiles string which encodes the molecule with the desired settings.
        """
        return RDKit.get_smiles(
            rdkit_mol=self.to_rdkit(),
            isomeric=isomeric,
            explicit_hydrogens=explicit_hydrogens,
            mapped=mapped,
        )

    def generate_atom_names(self) -> None:
        """
        Generate a unique set of atom names for the molecule.
        """
        atom_names = {}
        for atom in self.atoms:
            symbol = atom.atomic_symbol
            if symbol not in atom_names:
                atom_names[symbol] = 1
            else:
                atom_names[symbol] += 1

            atom.atom_name = f"{symbol}{atom_names[symbol]}"

    def _validate_atom_names(self) -> None:
        """
        Check that the ligand has unique atom names if not generate a new set.
        """
        if not self.has_unique_atom_names:
            self.generate_atom_names()

    def to_qcschema(
        self, extras: Optional[Dict] = None
    ) -> qcel.models.Molecule:
        """
        build a qcschema molecule from the ligand object, this is useful to interface with QCEngine and QCArchive.
        """
        import copy

        # make sure we have a conformer
        if self.coordinates == [] or self.coordinates is None:
            raise ConformerError(
                "The molecule must have a conformation to make a qcschema molecule."
            )
        coords = copy.deepcopy(self.coordinates)
        # input must be in bohr
        coords *= ANGS_TO_BOHR
        # we do not store explicit bond order so guess at 1
        bonds = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order)
            for bond in self.bonds
        ]
        mapped_smiles = self.to_smiles(
            isomeric=True, explicit_hydrogens=True, mapped=True
        )
        if extras is not None:
            extras["canonical_isomeric_explicit_hydrogen_mapped_smiles"] = (
                mapped_smiles
            )
        else:
            extras = {
                "canonical_isomeric_explicit_hydrogen_mapped_smiles": mapped_smiles
            }

        symbols = [atom.atomic_symbol for atom in self.atoms]
        schema_info = {
            "symbols": symbols,
            "geometry": coords,
            "connectivity": bonds,
            "molecular_charge": self.charge,
            "molecular_multiplicity": self.multiplicity,
            "extras": extras,
            "fix_com": True,
            "fix_orientation": True,
            "fix_symmetry": "c1",
        }
        return qcel.models.Molecule.from_data(schema_info, validate=True)

    def add_conformer(self, file_name: str) -> None:
        """
        Read the given input file extract  the conformers and save them to the ligand.
        TODO do we want to check that the connectivity is the same?
        """
        input_data = ReadInput.from_file(file_name=file_name)
        if input_data.coords is None:
            # get the coords from the rdkit molecule
            coords = input_data.rdkit_mol.GetConformer().GetPositions()
        else:
            if isinstance(input_data.coords, list):
                coords = input_data.coords[-1]
            else:
                coords = input_data.coords
        self.coordinates = coords


class ModSemMaths:
    """Static methods for various mathematical functions relevant to the modified Seminario method."""

    def __repr__(self):
        return f"{self.__class__.__name__}({self.__dict__!r})"

    @staticmethod
    def unit_vector_normal_to_bond(u_bc, u_ab):
        """Calculates unit vector which is normal to the plane abc."""

        cross = np.cross(u_bc, u_ab)

        return cross / np.linalg.norm(cross)

    @staticmethod
    def unit_vector_along_bond(coords, bond):
        """Calculates the unit vector along a bond."""

        atom_a, atom_b = bond
        diff_ab = coords[atom_b] - coords[atom_a]

        return diff_ab / np.linalg.norm(diff_ab)

    @staticmethod
    def u_pa_from_angles(angle, coords):
        """This gives the vector in the plane a, b, c and perpendicular to a to b."""

        atom_a, atom_b, atom_c = angle

        u_ab = ModSemMaths.unit_vector_along_bond(coords, (atom_a, atom_b))
        u_cb = ModSemMaths.unit_vector_along_bond(coords, (atom_c, atom_b))

        u_n = ModSemMaths.unit_vector_normal_to_bond(u_cb, u_ab)

        return ModSemMaths.unit_vector_normal_to_bond(u_n, u_ab)

    @staticmethod
    def dot_product(u_pa, eig_ab):

        return sum(u_pa[i] * eig_ab[i].conjugate() for i in range(3))

    @staticmethod
    def force_constant_bond(bond, eigenvals, eigenvecs, coords):
        """Force Constant - Equation 10 of Seminario paper - gives force constant for bond."""

        atom_a, atom_b = bond
        eigenvals_ab = eigenvals[atom_a, atom_b, :]
        eigenvecs_ab = eigenvecs[:, :, atom_a, atom_b]

        unit_vectors_ab = ModSemMaths.unit_vector_along_bond(coords, bond)

        return -0.5 * sum(
            eigenvals_ab[i] * abs(np.dot(unit_vectors_ab, eigenvecs_ab[:, i]))
            for i in range(3)
        )

    @staticmethod
    def force_constant_angle(
        angle, bond_lens, eigenvals, eigenvecs, coords, scalings
    ):
        """
        Force Constant - Equation 14 of Seminario paper - gives force constant for angle
        (in kcal/mol/rad^2) and equilibrium angle (in degrees).
        """

        atom_a, atom_b, atom_c = angle

        u_ab = ModSemMaths.unit_vector_along_bond(coords, (atom_a, atom_b))
        u_cb = ModSemMaths.unit_vector_along_bond(coords, (atom_c, atom_b))

        bond_len_ab = bond_lens[atom_a, atom_b]
        eigenvals_ab = eigenvals[atom_a, atom_b, :]
        eigenvecs_ab = eigenvecs[:3, :3, atom_a, atom_b]

        bond_len_bc = bond_lens[atom_b, atom_c]
        eigenvals_cb = eigenvals[atom_c, atom_b, :]
        eigenvecs_cb = eigenvecs[:3, :3, atom_c, atom_b]

        # Normal vector to angle plane found
        u_n = ModSemMaths.unit_vector_normal_to_bond(u_cb, u_ab)

        # Angle is linear:
        if abs(np.linalg.norm(u_cb - u_ab)) < 0.01 or (
            1.99 < abs(np.linalg.norm(u_cb - u_ab)) < 2.01
        ):
            # Scalings are set to 1.
            k_theta, theta_0 = ModSemMaths.f_c_a_special_case(
                u_ab,
                u_cb,
                [bond_len_ab, bond_len_bc],
                [eigenvals_ab, eigenvals_cb],
                [eigenvecs_ab, eigenvecs_cb],
            )

        else:
            u_pa = ModSemMaths.unit_vector_normal_to_bond(u_n, u_ab)
            u_pc = ModSemMaths.unit_vector_normal_to_bond(u_cb, u_n)

            # Scaling due to additional angles - Modified Seminario Part
            sum_first = (
                sum(
                    eigenvals_ab[i]
                    * abs(ModSemMaths.dot_product(u_pa, eigenvecs_ab[:, i]))
                    for i in range(3)
                )
                / scalings[0]
            )
            sum_second = (
                sum(
                    eigenvals_cb[i]
                    * abs(ModSemMaths.dot_product(u_pc, eigenvecs_cb[:, i]))
                    for i in range(3)
                )
                / scalings[1]
            )

            # Added as two springs in series
            k_theta = (1 / ((bond_len_ab**2) * sum_first)) + (
                1 / ((bond_len_bc**2) * sum_second)
            )
            k_theta = 1 / k_theta

            # Change to OPLS form
            k_theta = abs(k_theta * 0.5)

            # Equilibrium Angle
            theta_0 = np.degrees(np.arccos(np.dot(u_ab, u_cb)))

        return k_theta, theta_0

    @staticmethod
    def f_c_a_special_case(u_ab, u_cb, bond_lens, eigenvals, eigenvecs):
        """
        Force constant angle special case, for example nitrile groups.
        This is for when the bond is linear, and therefore cannot be sampled around in the same way.
        The perpendicular vector is not defined for a linear bond.
        """

        # Number of samples around the bond.
        n_samples = 200
        k_theta_array = np.zeros(n_samples)

        for theta in range(n_samples):

            u_n = [
                np.sin(theta) * np.cos(theta),
                np.sin(theta) * np.sin(theta),
                np.cos(theta),
            ]

            u_pa = ModSemMaths.unit_vector_normal_to_bond(u_n, u_ab)
            u_pc = ModSemMaths.unit_vector_normal_to_bond(u_cb, u_n)

            sum_first = sum(
                eigenvals[0][i]
                * abs(ModSemMaths.dot_product(u_pa, eigenvecs[0][:, i]))
                for i in range(3)
            )
            sum_second = sum(
                eigenvals[1][i]
                * abs(ModSemMaths.dot_product(u_pc, eigenvecs[1][:, i]))
                for i in range(3)
            )

            k_theta_i = (1 / ((bond_lens[0] ** 2) * sum_first)) + (
                1 / ((bond_lens[1] ** 2) * sum_second)
            )
            k_theta_i = 1 / k_theta_i

            k_theta_array[theta] = abs(k_theta_i * 0.5)

        k_theta = np.average(k_theta_array)
        theta_0 = np.degrees(np.arccos(np.dot(u_ab, u_cb)))

        return k_theta, theta_0


class ModSeminario(StageBase):

    type: Literal["ModSeminario"] = "ModSeminario"
    vibrational_scaling: float = Field(
        1.0,
        description="The vibration scaling that should be used to correct the reference DFT frequencies.",
    )

    @classmethod
    def is_available(cls) -> bool:
        """This class is part of qubekit and always available."""
        return True

    def start_message(self, **kwargs) -> str:
        return "Calculating new bond and angle parameters with the modified Seminario method."

    def finish_message(self, **kwargs) -> str:
        return "Bond and angle parameters calculated."

    def run(self, molecule: Ligand, **kwargs) -> Ligand:
        """
        The main worker stage which takes the molecule and its hessian and calculates the modified seminario method.

        Args:
            molecule: The qubekit molecule class that should contain a valid hessian and optimised coordinates.

        Note:
            Please cite this method using <J. Chem. Theory Comput. (2018), doi:10.1021/acs.jctc.7b00785>
        """

        # reset the bond and angle parameter groups
        molecule.BondForce.clear_parameters()
        molecule.AngleForce.clear_parameters()
        # convert the hessian from atomic units
        conversion = HA_TO_KCAL_P_MOL / (BOHR_TO_ANGS**2)
        # make sure we do not change the molecule hessian
        hessian = copy.deepcopy(molecule.hessian)
        hessian *= conversion
        self._modified_seminario_method(molecule=molecule, hessian=hessian)
        # apply symmetry to the bond and angle parameters
        molecule.symmetrise_bonded_parameters()

        return molecule

    def _modified_seminario_method(
        self, molecule: Ligand, hessian: np.ndarray
    ) -> Ligand:
        """
        Calculate the new bond and angle terms after being passed the symmetric Hessian and
        optimised molecule coordinates.
        """
        size_mol = molecule.n_atoms
        eigenvecs = np.empty((3, 3, size_mol, size_mol), dtype=complex)
        eigenvals = np.empty((size_mol, size_mol, 3), dtype=complex)
        bond_lens = np.zeros((size_mol, size_mol))

        for i in range(size_mol):
            for j in range(size_mol):
                diff_i_j = (
                    molecule.coordinates[i, :] - molecule.coordinates[j, :]
                )
                bond_lens[i, j] = np.linalg.norm(diff_i_j)

                partial_hessian = hessian[
                    (i * 3) : ((i + 1) * 3), (j * 3) : ((j + 1) * 3)
                ]

                eigenvals[i, j, :], eigenvecs[:, :, i, j] = np.linalg.eig(
                    partial_hessian
                )

        # The bond and angle values are calculated and written to file.
        self.calculate_bonds(eigenvals, eigenvecs, molecule, bond_lens)
        self.calculate_angles(eigenvals, eigenvecs, molecule, bond_lens)
        return molecule

    def calculate_angles(
        self, eigenvals, eigenvecs, molecule: Ligand, bond_lengths: np.ndarray
    ):
        """
        Uses the modified Seminario method to find the angle parameters and prints them to file.
        """

        # A structure is created with the index giving the central atom of the angle;
        # an array then lists the angles with that central atom.
        # e.g. central_atoms_angles[3] contains an array of angles with central atom 3.

        # Connectivity information for Modified Seminario Method
        central_atoms_angles = []

        for coord in range(molecule.n_atoms):
            central_atoms_angles.append([])
            for count, angle in enumerate(molecule.angles):
                if coord == angle[1]:
                    # For angle abc, atoms a, c are written to array
                    central_atoms_angles[coord].append(
                        [angle[0], angle[2], count]
                    )

                    # For angle abc, atoms c a are written to array
                    central_atoms_angles[coord].append(
                        [angle[2], angle[0], count]
                    )

        # Sort rows by atom number
        for coord in range(molecule.n_atoms):
            central_atoms_angles[coord] = sorted(
                central_atoms_angles[coord], key=itemgetter(0)
            )

        # Find normals u_pa for each angle
        unit_pa_all_angles = []

        for i in range(len(central_atoms_angles)):
            unit_pa_all_angles.append([])
            for j in range(len(central_atoms_angles[i])):
                # For the angle at central_atoms_angles[i][j,:] the u_pa value is found for plane abc and bond ab,
                # where abc corresponds to the order of the arguments. This is why the reverse order was also added.
                angle = (
                    central_atoms_angles[i][j][0],
                    i,
                    central_atoms_angles[i][j][1],
                )
                unit_pa_all_angles[i].append(
                    ModSemMaths.u_pa_from_angles(angle, molecule.coordinates)
                )

        # Finds the contributing factors from the other angle terms
        scaling_factor_all_angles = []

        for i in range(len(central_atoms_angles)):
            scaling_factor_all_angles.append([])
            for j in range(len(central_atoms_angles[i])):
                n = m = 1
                angles_around = extra_contribs = 0
                scaling_factor_all_angles[i].append([0, 0])

                # Position in angle list
                scaling_factor_all_angles[i][j][1] = central_atoms_angles[i][
                    j
                ][2]

                # Goes through the list of angles with the same central atom, then computes the term needed for MSM.

                # Forwards direction, finds the same bonds with the central atom i
                while (
                    (j + n) < len(central_atoms_angles[i])
                ) and central_atoms_angles[i][j][0] == central_atoms_angles[i][
                    j + n
                ][
                    0
                ]:
                    extra_contribs += (
                        abs(
                            np.dot(
                                unit_pa_all_angles[i][j][:],
                                unit_pa_all_angles[i][j + n][:],
                            )
                        )
                    ) ** 2
                    n += 1
                    angles_around += 1

                # Backwards direction, finds the same bonds with the central atom i
                while ((j - m) >= 0) and central_atoms_angles[i][j][
                    0
                ] == central_atoms_angles[i][j - m][0]:
                    extra_contribs += (
                        abs(
                            np.dot(
                                unit_pa_all_angles[i][j][:],
                                unit_pa_all_angles[i][j - m][:],
                            )
                        )
                    ) ** 2
                    m += 1
                    angles_around += 1

                scaling_factor_all_angles[i][j][0] = 1
                if n != 1 or m != 1:
                    # Finds the mean value of the additional contribution
                    scaling_factor_all_angles[i][j][0] += extra_contribs / (
                        m + n - 2
                    )

        scaling_factors_angles_list = [[] for _ in range(molecule.n_angles)]

        # Orders the scaling factors according to the angle list
        for i in range(len(central_atoms_angles)):
            for j in range(len(central_atoms_angles[i])):
                scaling_factors_angles_list[
                    scaling_factor_all_angles[i][j][1]
                ].append(scaling_factor_all_angles[i][j][0])

        k_theta, theta_0 = np.zeros(len(molecule.angles)), np.zeros(
            len(molecule.angles)
        )

        conversion = KCAL_TO_KJ * 2

        with open("Modified_Seminario_Angles.txt", "w") as angle_file:

            for i, angle in enumerate(molecule.angles):

                scalings = scaling_factors_angles_list[i]

                # Ensures that there is no difference when the ordering is changed.
                ab_k_theta, ab_theta_0 = ModSemMaths.force_constant_angle(
                    angle,
                    bond_lengths,
                    eigenvals,
                    eigenvecs,
                    molecule.coordinates,
                    scalings,
                )
                ba_k_theta, ba_theta_0 = ModSemMaths.force_constant_angle(
                    angle[::-1],
                    bond_lengths,
                    eigenvals,
                    eigenvecs,
                    molecule.coordinates,
                    scalings[::-1],
                )

                # Vib_scaling takes into account DFT deficiencies / anharmonicity.
                k_theta[i] = ((ab_k_theta + ba_k_theta) / 2) * (
                    self.vibrational_scaling**2
                )
                theta_0[i] = (ab_theta_0 + ba_theta_0) / 2

                angle_file.write(
                    f"{molecule.atoms[angle[0]].atom_name}-{molecule.atoms[angle[1]].atom_name}-{molecule.atoms[angle[2]].atom_name}  "
                )
                angle_file.write(
                    f"{k_theta[i]:.3f}   {theta_0[i]:.3f}   {angle[0]}   {angle[1]}   {angle[2]}\n"
                )

                # Add ModSem values to ligand object.
                molecule.AngleForce.create_parameter(
                    atoms=angle,
                    angle=theta_0[i] * DEG_TO_RAD,
                    k=k_theta[i] * conversion,
                )

    def calculate_bonds(
        self, eigenvals, eigenvecs, molecule: Ligand, bond_lengths: np.ndarray
    ):
        """
        Uses the modified Seminario method to find the bond parameters and print them to file.
        """

        bonds = molecule.to_topology().edges
        conversion = KCAL_TO_KJ * 200

        k_b, bond_len_list = np.zeros(len(bonds)), np.zeros(len(bonds))

        with open("Modified_Seminario_Bonds.txt", "w") as bond_file:

            for pos, bond in enumerate(bonds):
                ab = ModSemMaths.force_constant_bond(
                    bond, eigenvals, eigenvecs, molecule.coordinates
                )
                ba = ModSemMaths.force_constant_bond(
                    bond[::-1], eigenvals, eigenvecs, molecule.coordinates
                )

                # Order of bonds sometimes causes slight differences; find the mean and apply vib_scaling.
                k_b[pos] = np.real((ab + ba) / 2) * (
                    self.vibrational_scaling**2
                )

                bond_len_list[pos] = bond_lengths[bond]
                bond_file.write(
                    f"{molecule.atoms[bond[0]].atom_name}-{molecule.atoms[bond[1]].atom_name}  "
                )
                bond_file.write(
                    f"{k_b[pos]:.3f}   {bond_len_list[pos]:.3f}   {bond[0]}   {bond[1]}\n"
                )

                # Add ModSem values to ligand object.
                molecule.BondForce.create_parameter(
                    atoms=bond,
                    length=bond_len_list[pos] / 10,
                    k=conversion * k_b[pos],
                )
