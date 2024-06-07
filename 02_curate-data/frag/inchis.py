# generate a list of inchikeys for all of the molecules in our current training
# and testing data sets to filter out when querying the database

import logging
from pathlib import Path

from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit import Molecule
from tqdm import tqdm

logging.getLogger("openff").setLevel(logging.ERROR)

vf = Path("/home/brent/omsf/projects/valence-fitting")

datasets = [
    vf / "02_curate-data/datasets/combined-opt.json",
    vf / "02_curate-data/datasets/combined-td.json",
    "/home/brent/omsf/projects/benchmarking/datasets/industry.json",
]

inchis = set()
for ds in datasets:
    if "td" in str(ds):
        d = TorsionDriveResultCollection.parse_file(ds)
    else:
        d = OptimizationResultCollection.parse_file(ds)
    values = d.entries.values()
    assert len(values) == 1
    records = list(values)[0]
    for cmiles in tqdm((rec.cmiles for rec in records), total=len(records)):
        inchis.add(
            Molecule.from_smiles(
                cmiles, allow_undefined_stereo=True
            ).to_inchikey()
        )

outfile = "inchis.dat"
print(f"writing {len(inchis)} inchis to {outfile}")
with open(outfile, "w") as out:
    for inchi in inchis:
        print(inchi, file=out)
