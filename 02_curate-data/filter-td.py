import json
from curate_dataset import filter_td_data
from openff.qcsubmit.results import TorsionDriveResultCollection

dataset = TorsionDriveResultCollection(entries=dict(json.loads(r"""{json}""")))
dataset = filter_td_data(
    dataset,
    "td_records_to_remove.dat",
    include_iodine=False,
)
print(dataset.json())
