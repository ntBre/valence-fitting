from dataclasses import dataclass
from typing import Iterator

import click
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit import ForceField, Molecule
from tqdm import tqdm

from main import param_sort_key


@dataclass
class Match:
    env: int
    rec: set[str]
    mol: set[str]
    tor: int

    def __init__(self, env=0, rec=None, mol=None, tor=0):
        self.env = env
        self.tor = tor
        if rec is None:
            self.rec = set()
        if mol is None:
            self.mol = set()


def load_dataset(ds) -> Iterator[tuple[str, str]]:
    "Yields a sequence of record id, smiles pairs"
    seen = set()
    for rec in next(iter(ds.entries.values())):
        if (s := rec.cmiles) not in seen:
            seen.add(s)
            yield rec.record_id, s


def td_main(dataset, ff):
    ds = TorsionDriveResultCollection.parse_file(dataset)
    print("calling to_records")
    for rec, mol in tqdm(ds.to_records()):
        print(rec.record_id, mol.to_smiles())


def opt_main(dataset, ff):
    ds = load_dataset(OptimizationResultCollection.parse_file(dataset))
    # pre-fill instead of defaultdict to catch zeros
    matches = {
        p.id: Match()
        for p in ff.get_parameter_handler("ProperTorsions").parameters
    }
    for rec_id, smiles in ds:
        mol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        labels = ff.label_molecules(mol.to_topology())[0]["ProperTorsions"]
        for _env, p in labels.items():
            matches[p.id].env += 1
            matches[p.id].rec.add(rec_id)
            matches[p.id].mol.add(smiles)

    matches = list(sorted(matches.items(), key=lambda x: param_sort_key(x[0])))
    pid_w, env_w, rec_w, smi_w = 6, 8, 8, 8
    print(
        f"{'pid':<{pid_w}} {'env':>{env_w}} {'rec':>{rec_w}} {'smi':>{smi_w}}"
    )
    for pid, m in matches:
        rec, smi = len(m.rec), len(m.mol)
        print(f"{pid:<{pid_w}} {m.env:>{env_w}} {rec:>{rec_w}} {smi:>{smi_w}}")


@click.command()
@click.option("--forcefield", "-f")
@click.option("--dataset", "-d")
@click.option("--torsions", "-t", is_flag=True)
def main(forcefield, dataset, torsions):
    ff = ForceField(forcefield)
    if torsions:
        td_main(dataset, ff)
    else:
        opt_main(dataset, ff)
    pass


if __name__ == "__main__":
    main()
