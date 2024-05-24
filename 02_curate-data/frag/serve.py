import re

from flask import Flask
from jinja2 import Environment, PackageLoader, select_autoescape

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
