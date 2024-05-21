import json

from query import Match

with open("out.query.json") as f:
    data: dict[str, dict] = json.load(f)

matches: dict[str, Match] = {k: Match(**v) for k, v in data.items()}

for m in matches.values():
    lm = len(m.molecules)
    print(m.pid, m.smirks, lm)
    if lm <= 5:
        for mol in m.molecules:
            print("\t", mol)
