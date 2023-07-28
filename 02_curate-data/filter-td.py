import json
from curate_dataset import filter_td_data
from openff.qcsubmit.results import TorsionDriveResultCollection

d = json.loads(r"""{json}""")
dataset = TorsionDriveResultCollection(entries=d['entries'])
dataset = filter_td_data(
    dataset,
    "td_records_to_remove.dat",
    include_iodine=False,
)
print(dataset.json())
