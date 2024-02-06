# prints the difference in torsion parameters between the TM force field and
# the original Sage 2.1.0 force field

import re
from sys import argv

from openff.toolkit import ForceField

# pasted from benchmarking/parse_hist.py
LABEL = re.compile(r"([bati])(\d+)([a-z]*)")


def sort_label(key):
    t, n, tail = LABEL.match(key).groups()
    return (t, int(n), tail)


def main():
    if len(argv) < 3:
        print("Usage: diff ff1 ff2")
        exit(1)
    ff1 = ForceField(argv[1], allow_cosmetic_attributes=True)
    ff2 = ForceField(argv[2], allow_cosmetic_attributes=True)

    h1 = ff1.get_parameter_handler("Bonds")
    h2 = ff2.get_parameter_handler("Bonds")
    print(f"{'Pid':5} {'k1':>8} {'k2':>8} {'L1':>5} {'L2':>5}")
    for t1 in h1.parameters:
        t2 = h2[t1.smirks]
        b1k, b2k = t1.k.magnitude, t2.k.magnitude
        b1l, b2l = t1.length.magnitude, t2.length.magnitude
        print(f"{t1.id:5} {b1k:8.2f} {b2k:8.2f} {b1l:5.2f} {b2l:5.2f}")

    print()

    h1 = ff1.get_parameter_handler("Angles")
    h2 = ff2.get_parameter_handler("Angles")
    print(f"{'Pid':5} {'k1':>8} {'k2':>8} {'θ1':>8} {'θ2':>8}")
    for t1 in h1.parameters:
        t2 = h2[t1.smirks]
        b1k, b2k = t1.k.magnitude, t2.k.magnitude
        b1l, b2l = t1.angle.magnitude, t2.angle.magnitude
        print(f"{t1.id:5} {b1k:8.2f} {b2k:8.2f} {b1l:8.2f} {b2l:8.2f}")

    print()

    h1 = ff1.get_parameter_handler("ProperTorsions")
    h2 = ff2.get_parameter_handler("ProperTorsions")

    print(f"{'Pid':5} ", end="")
    for i in range(6):
        print(f"{f'k{i+1}-1':>8} {f'k{i+1}-2':>8}", end=" ")
    print()
    for t1 in h1.parameters:
        t2 = h2[t1.smirks]
        print(f"{t1.id:5}", end=" ")
        for i in range(6):
            k = f"k{i+1}"
            b1k = getattr(t1, k, None)
            if b1k is None:
                continue
            b1k, b2k = b1k.magnitude, getattr(t2, k).magnitude
            print(f"{b1k:8.2f} {b2k:8.2f}", end=" ")
        print()


if __name__ == "__main__":
    main()
