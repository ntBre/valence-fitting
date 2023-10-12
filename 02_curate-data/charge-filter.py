import click
from vflib import Molecules, load_dataset

from filters import ChargeCheckFilter


@click.command()
@click.option("--input")
@click.option("--output")
def main(input, output):
    dataset = load_dataset(input)
    init = len(Molecules(dataset))
    dataset = dataset.filter(ChargeCheckFilter())
    with open(output, "w") as out:
        out.write(dataset.json(indent=2))
    final = len(Molecules(dataset))
    print(f"filtered {init} entries in {input} to {final}")


if __name__ == "__main__":
    main()
