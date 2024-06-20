from collections import defaultdict

d = defaultdict(list)
with open("opt.smi") as f:
    for line in f:
        if line.startswith("#"):
            continue
        pid, smiles, *_ = line.split()
        d[pid].append(smiles)


with open("train.opt.smi", "w") as train, open("bench.opt.smi", "w") as bench:
    for k, v in d.items():
        for smiles in v[:-1]:  # put all but one in train
            print(k, smiles, file=train)
        print(k, v[-1], file=bench)


with open("train.opt.smi") as train, open("bench.opt.smi") as bench:
    lines = [line for line in train]
    lines.extend([line for line in bench])

assert len(lines) == sum((len(v) for v in d.values()))
