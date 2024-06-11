with open("want.params") as inp:
    want = {line.strip() for line in inp}

with open("dataset.smi") as inp:
    got = {line.strip().split()[0] for line in inp}

with open("remaining.dat", "w") as out:
    for p in want - got:
        print(p, file=out)
