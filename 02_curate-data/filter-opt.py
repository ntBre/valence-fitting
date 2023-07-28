import json
from curate_dataset import filter_opt_data
from openff.qcsubmit.results import OptimizationResultCollection

d = json.loads(r"""{json}""")
dataset = OptimizationResultCollection(entries=d['entries'])
dataset = filter_opt_data(
    dataset,
    "opt_records_to_remove.dat",
    include_iodine=False,
)
print(dataset.json())
