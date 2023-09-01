# first check how much data we include from Sage vs new data

import json


def load(filename):
    with open(filename, "r") as inp:
        res = json.load(inp)
    return [v for vals in res["entries"].values() for v in vals]


def to_set(d):
    return {r["record_id"] for r in d}


sage_opt = load("02_curate-data/sage/opt.json")
sage_td = load("02_curate-data/sage/td.json")
pavan_opt = load("02_curate-data/datasets/core-opt.json")
pavan_td = load("02_curate-data/datasets/core-td.json")

sage_opt_records = to_set(sage_opt)
sage_td_records = to_set(sage_td)
pavan_opt_records = to_set(pavan_opt)
pavan_td_records = to_set(pavan_td)

print()
print("Before Processing")
print("-----------------")
print(f"{len(sage_opt_records)} Opt records from Sage")
print(f"{len(sage_td_records)} TD records from Sage")
print(f"{len(pavan_opt_records)} Opt records from Pavan")
print(f"{len(pavan_td_records)} TD records from Pavan")

combined_opt_records = sage_opt_records | pavan_opt_records
combined_td_records = sage_td_records | pavan_td_records

print()
print("After Combining")
print("---------------")
print(f"{len(combined_opt_records)} Opt records combined")
print(f"{len(combined_td_records)} TD records combined")

assert len(combined_opt_records) == len(sage_opt_records) + len(
    pavan_opt_records
)
assert len(combined_td_records) == len(sage_td_records) + len(pavan_td_records)

# we can also double-check that the combined JSON files from the
# valence-fitting pipeline have these numbers as well
my_combined_opt = load("02_curate-data/datasets/combined-opt.json")
my_combined_td = load("02_curate-data/datasets/combined-td.json")
my_combined_opt = to_set(my_combined_opt)
my_combined_td = to_set(my_combined_td)

print()
print("Combining Sanity Check")
print("----------------------")
print(f"{len(my_combined_opt)} Opt records combined")
print(f"{len(my_combined_td)} TD records combined")

assert len(my_combined_opt) == len(combined_opt_records)
assert len(my_combined_td) == len(combined_td_records)

# so all of the data in the TM set is distinct from that in Sage. now we need
# to check what gets filtered out

filtered_opt = to_set(load("02_curate-data/datasets/filtered-opt.json"))
filtered_td = to_set(load("02_curate-data/datasets/filtered-td.json"))

# how many filtered records were left, and how many of them were from Sage vs
# from TM?
print()
print("Filtered")
print("--------")
print(f"{len(filtered_opt)} Opt records")
print(f"\t{len(filtered_opt & sage_opt_records)} from Sage")
print(f"\t{len(filtered_opt & pavan_opt_records)} from TM")
print(f"{len(filtered_td)} TD records")
print(f"\t{len(filtered_td & sage_td_records)} from Sage")
print(f"\t{len(filtered_td & pavan_td_records)} from TM")

# Nearly half of the TM data for both opt and td is filtered out. On the one
# hand, this raises the question of why bother including it, but on the other
# hand, it makes it surprising that this small amount of data worsens the force
# field.
