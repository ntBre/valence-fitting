import re
from dataclasses import dataclass

from flask import Flask, request, send_from_directory
from jinja2 import Environment, PackageLoader, select_autoescape
from openff.toolkit import ForceField
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage, rdDepictor, rdMolDraw2D

from query import find_matches, into_params, mol_from_smiles
from store import Store

app = Flask("serve")
env = Environment(
    loader=PackageLoader("serve"), autoescape=select_autoescape()
)
dbname = "store.sqlite"
ffname = (
    "../../01_generate-forcefield/output/"
    "initial-force-field-openff-2.1.0.offxml"
)
off = ForceField(ffname)
pid_to_smarts = {
    p.id: p.smirks
    for p in off.get_parameter_handler("ProperTorsions").parameters
}


PID_RE = re.compile("(t)([0-9]+)(.*)")


def pid_sort(pid: str) -> tuple[str, int, str | None]:
    "Return the fields of a ProperTorsion parameter ID as a tuple for sorting"
    t, n, tail = PID_RE.match(pid).groups()
    return (t, int(n), tail)


def mol_to_svg(
    mol: Chem.Mol, width: int, height: int, legend: str, hl_atoms: list[int]
) -> str:
    "Return an SVG of `mol`"
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    rdmol = rdMolDraw2D.PrepareMolForDrawing(mol)
    return MolsToGridImage(
        [rdmol],
        useSVG=True,
        highlightAtomLists=[hl_atoms],
        subImgSize=(300, 300),
        molsPerRow=1,
    )


@app.route("/js/<path:name>")
def js(name):
    return send_from_directory("js", name)


@app.route("/css/<path:name>")
def css(name):
    return send_from_directory("css", name)


@app.route("/")
def index():
    template = env.get_template("index.html")
    table = Store.quick()
    pairs = [
        (m.pid, len(m.molecules)) for m in table.get_forcefield(ffname).matches
    ]
    pairs = sorted(pairs, key=lambda x: pid_sort(x[0]))
    parameter_ids, molecule_counts = (
        [x[0] for x in pairs],
        [x[1] for x in pairs],
    )
    return template.render(
        parameter_ids=parameter_ids,
        molecule_counts=molecule_counts,
        pid_counts=[0] * len(molecule_counts),  # TODO get from db
    )


@dataclass
class DrawMol:
    smiles: str
    natoms: int
    svg: str


@app.route("/param/<pid>")
def param(pid):
    MAX_DRAW = 50
    max_draw = int(request.args.get("max", MAX_DRAW))
    table = Store.quick()
    mols = get_smiles_list(table, ffname, pid)
    mol_map = into_params(off)

    draw_mols = mols_to_draw(mols, pid, mol_map, max_draw)

    smarts = pid_to_smarts[pid]  # TODO error page on missing
    total_mols = len(mols)
    template = env.get_template("param.html")
    return template.render(
        smarts=smarts,
        pid=pid,
        mols=draw_mols,
        cur_mols=len(draw_mols),
        total_mols=total_mols,
    )


def get_smiles_list(table, ffname, pid) -> list[tuple[Chem.Mol, str, int]]:
    """Return a sequence of Mol, SMILES, natom tuples from `table`, sorted by
    natoms.

    """
    smiles_list = table.get_smiles_matching(ffname, pid)
    mols = []
    for s in smiles_list:
        mol = mol_from_smiles(s)
        natoms = mol.GetNumAtoms()
        mols.append((mol, s, natoms))
    mols = sorted(mols, key=lambda k: k[2])  # sort by natoms
    return mols


def mols_to_draw(mols, pid, mol_map, max_draw):
    draw_mols = []
    for mol, smiles, natoms in mols[:max_draw]:
        matches = find_matches(mol_map, mol)
        hl_atoms = []
        for atoms, mpid in matches.items():
            if mpid == pid:
                hl_atoms = atoms
                break
        svg = mol_to_svg(mol, 300, 300, "", hl_atoms)
        draw_mols.append(DrawMol(smiles, natoms, svg))
    return draw_mols
