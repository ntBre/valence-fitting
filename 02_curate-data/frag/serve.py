import re

from flask import Flask, request
from jinja2 import Environment, PackageLoader, select_autoescape
from openff.toolkit import ForceField

from query import mol_from_smiles
from store import Store

app = Flask("serve")
env = Environment(
    loader=PackageLoader("serve"), autoescape=select_autoescape()
)
dbname = "store.sqlite"
ff = (
    "../../01_generate-forcefield/output/"
    "initial-force-field-openff-2.1.0.offxml"
)
off = ForceField(ff)
pid_to_smarts = {
    p.id: p.smirks
    for p in off.get_parameter_handler("ProperTorsions").parameters
}


PID_RE = re.compile("(t)([0-9]+)(.*)")


def pid_sort(pid: str) -> tuple[str, int, str | None]:
    "Return the fields of a ProperTorsion parameter ID as a tuple for sorting"
    t, n, tail = PID_RE.match(pid).groups()
    return (t, int(n), tail)


@app.route("/")
def index():
    template = env.get_template("index.html")
    table = Store.quick()
    pairs = [
        (m.pid, len(m.molecules)) for m in table.get_forcefield(ff).matches
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


@app.route("/param/<pid>")
def param(pid):
    MAX_DRAW = 50
    max_draw = request.args.get("max", MAX_DRAW)
    table = Store.quick()
    smiles_list = table.get_smiles_matching(ff, pid)
    mols = []
    for s in smiles_list:
        mol = mol_from_smiles(s)
        natoms = mol.GetNumAtoms()
        mols.append((mol, s, natoms))
    mols = sorted(mols, key=lambda k: k[2])  # sort by natoms
    smarts = pid_to_smarts[pid]  # TODO error page on missing
    total_mols = len(mols)
    template = env.get_template("param.html")
    return template.render(
        smarts=smarts,
        pid=pid,
        mols=mols,
        cur_mols=len(mols),
        total_mols=total_mols,
    )
