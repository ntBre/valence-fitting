import re
import warnings
from dataclasses import dataclass
from http import HTTPStatus

from flask import Flask, request, send_from_directory
from jinja2 import Environment, PackageLoader, select_autoescape
from openff.toolkit import ForceField
from rdkit import Chem

from parse import tanimoto
from query import find_matches, into_params, mol_from_smiles
from store import Store
from utils import find_smallest, make_svg, mol_to_svg

warnings.filterwarnings("ignore")
with warnings.catch_warnings():
    from sklearn import cluster as skcluster

app = Flask("serve")
env = Environment(
    loader=PackageLoader("serve"), autoescape=select_autoescape()
)
env.globals["Chem"] = Chem
env.globals["make_svg"] = make_svg
env.globals["find_smallest"] = find_smallest
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


@dataclass
class Report:
    max: int
    nfps: int
    noise: int
    clusters: list[list[int]]
    mols: list[Chem.Mol]
    map: dict[str, str]
    mol_map: list[tuple[str, Chem.Mol]]


def make_fps(mols: list[Chem.Mol], radius: int):
    fpgen = Chem.AllChem.GetMorganGenerator(radius=radius)
    return [fpgen.GetFingerprint(mol) for mol in mols]


def dbscan(dist, eps, min_pts):
    db = skcluster.DBSCAN(eps=eps, min_samples=min_pts)
    return db.fit(dist).labels_


def make_cluster_report(ffname, mols, eps, min_pts) -> Report:
    MORGAN_RADIUS = 4
    map = {
        p.id: p.smirks
        for p in off.get_parameter_handler("ProperTorsions").parameters
    }
    mol_map = into_params(off)
    fps = make_fps(mols, MORGAN_RADIUS)
    nfps = len(fps)
    dist = tanimoto(fps)
    labels = dbscan(dist, eps, min_pts)
    label_max = max(labels)
    clusters = [[] for _ in range(label_max + 1)]
    noise = []
    for i, el in enumerate(labels):
        if el < 0:
            noise.append(i)
        else:
            clusters[el].append(i)

    noise_pts = len(noise)
    if len(clusters[0]) == 0:
        print(f"warning: all noise points: {labels}")
        clusters[0] = noise

    return Report(
        max=label_max,
        nfps=nfps,
        noise=noise_pts,
        clusters=clusters,
        mols=mols,
        map=map,
        mol_map=mol_map,
    )


@app.route("/cluster/<pid>")
def cluster(pid):
    smarts = pid_to_smarts[pid]
    table = Store.quick()
    mols = [mol_from_smiles(s) for s in table.get_smiles_matching(ffname, pid)]
    DBSCAN_EPS = 0.5
    DBSCAN_MIN_PTS = 1

    eps = request.args.get("eps", default=DBSCAN_EPS, type=float)
    min_pts = request.args.get("min_pts", default=DBSCAN_MIN_PTS, type=int)

    r = make_cluster_report(ffname, mols, eps, min_pts)

    clusters = sorted(r.clusters, key=lambda c: r.mols[c[0]].GetNumAtoms())
    template = env.get_template("cluster.html")
    return template.render(
        pid=pid,
        smarts=smarts,
        eps=eps,
        min_pts=min_pts,
        max=r.max,
        nfps=r.nfps,
        noise=r.noise,
        clusters=clusters,
        mols=r.mols,
        map=r.map,
        mol_map=r.mol_map,
    )


@app.route("/add-molecule", methods=["POST"])
def add_molecule():
    data = request.get_json()
    table = Store.quick()
    table.add_to_dataset(data["smiles"], data["pid"])
    return "", HTTPStatus.CREATED


@app.route("/preview-dataset")
def preview_dataset():
    return "preview dataset"