from openff.toolkit import ForceField
from openff.qcsubmit.results import TorsionDriveResultCollection
from collections import defaultdict
from sys import argv
import time


def check_coverage(
    ff="03_generate-initial-ff/output/initial-force-field-msm.offxml",
):
    start = time.time()

    ff = ForceField(
        ff,
        allow_cosmetic_attributes=True,
    )
    td_data = TorsionDriveResultCollection.parse_file(
        "02_curate-data/output/combined-td.json"
    )

    print(f"finished loading collection after {time.time() - start:.1f} sec")

    h = ff.get_parameter_handler("ProperTorsions")
    tors_ids = [p.id for p in h.parameters]

    records_and_molecules = td_data.to_records()

    print(f"finished to_records after {time.time() - start:.1f} sec")

    results = defaultdict(int)
    for _, molecule in records_and_molecules:
        all_labels = ff.label_molecules(molecule.to_topology())[0]
        torsions = all_labels["ProperTorsions"]
        for torsion in torsions.values():
            results[torsion.id] += 1

    print(f"finished counting results after {time.time() - start:.1f} sec")

    got = len(results)
    want = len(tors_ids)
    pct = 100.0 * float(got) / float(want)
    print(f"{got} / {want} ({pct:.1f}%) tors ids covered:")

    for id in tors_ids:
        smirk = h.get_parameter(dict(id=id))[0].smirks
        print(f"{id:5}{results[id]:5}   {smirk}")

    missing_ids = [k for k in results.keys() if results[k] == 0]
    missing_smirks = [
        h.get_parameter(dict(id=p))[0].smirks for p in missing_ids
    ]
    print("\nmissing ids:")
    for i, (id, smirk) in enumerate(zip(missing_ids, missing_smirks)):
        print(f"{i:5}{id:>7}   {smirk}")

    print(f"finished after {time.time() - start:.1f} sec")


if __name__ == "__main__":
    if len(argv) > 1:
        check_coverage(argv[1])
    else:
        check_coverage()
