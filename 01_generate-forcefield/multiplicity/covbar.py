import click
import matplotlib.pyplot as plt
import numpy as np


@click.command()
@click.option("--inp", "-i")
@click.option("--title", "-t")
@click.option("--output", "-o")
def main(inp, title, output):
    with open(inp) as inp:
        inp.readline()  # skip header
        labels, res = [], []
        for line in inp:
            [pid, env, rec, smi, tor] = line.split()
            labels.append(pid)
            res.append(int(tor))

    mx = np.max(res)
    bins = np.linspace(0, mx, num=mx + 1)
    plt.hist(res, bins=bins)
    plt.title(title)
    plt.xlabel("Number of matching torsion drives")
    plt.ylabel("Count")
    plt.savefig(output)


if __name__ == "__main__":
    main()
