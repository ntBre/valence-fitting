from openff.toolkit import ForceField
from openff.qcsubmit.results import TorsionDriveResultCollection
from collections import defaultdict

ff = ForceField(
    "03_generate-initial-ff/output/initial-force-field-msm.offxml",
    allow_cosmetic_attributes=True,
)
td_data = TorsionDriveResultCollection.parse_file(
    "02_curate-data/output/combined-td.json"
)

h = ff.get_parameter_handler("ProperTorsions")
tors_ids = [p.id for p in h.parameters]

records_and_molecules = td_data.to_records()

results = defaultdict(int)
for _, molecule in records_and_molecules:
    all_labels = ff.label_molecules(molecule.to_topology())[0]
    torsions = all_labels["ProperTorsions"]
    for torsion in torsions.values():
        results[torsion.id] += 1

got = len(results)
want = len(tors_ids)
print(
    f"{got} / {want} ({100.0 * float(got)/float(want):.1f}%) tors ids covered:"
)

for id in tors_ids:
    print(f"{id:5}{results[id]:5}")

missing_ids = [k for k in results.keys() if results[k] == 0]
missing_smirks = [h.get_parameter(dict(id=p))[0].smirks for p in missing_ids]
print("\nmissing ids:")
for i, (id, smirk) in enumerate(zip(missing_ids, missing_smirks)):
    print(f"{i:5}{id:>7}   {smirk}")
