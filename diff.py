# prints the difference in torsion parameters between the TM force field and
# the original Sage 2.1.0 force field

import re

from openff.toolkit import ForceField

# pasted from benchmarking/parse_hist.py
LABEL = re.compile(r"([bati])(\d+)([a-z]*)")


def sort_label(key):
    t, n, tail = LABEL.match(key).groups()
    return (t, int(n), tail)


s = ForceField("openff-2.1.0.offxml", allow_cosmetic_attributes=True)
t = ForceField(
    "01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml",
    allow_cosmetic_attributes=True,
)

s = s.get_parameter_handler("ProperTorsions")
t = t.get_parameter_handler("ProperTorsions")

s = set([p.id for p in s.parameters])
t = set([p.id for p in t.parameters])

for d in sorted(s.symmetric_difference(t), key=sort_label):
    print(d)
