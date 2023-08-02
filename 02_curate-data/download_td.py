from qcportal import FractalClient
from openff.qcsubmit.results import TorsionDriveResultCollection

ds_cache = "datasets/core-td.json"
td_datasets = ["OpenFF multiplicity correction torsion drive data v1.1"]

client = FractalClient()
dataset = TorsionDriveResultCollection.from_server(
    client=client,
    datasets=td_datasets,
    spec_name="default",
)
with open(ds_cache, "w") as out:
    out.write(dataset.json(indent=2))
