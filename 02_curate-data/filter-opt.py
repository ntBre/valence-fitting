import click
from openff.qcsubmit.results import OptimizationResultCollection

from filters import filter_opt_data


@click.command()
@click.option("-i/--input")
@click.option("-o/--output")
def main(input, output):
    dataset = OptimizationResultCollection.parse_file(input)
    dataset = filter_opt_data(
        dataset,
        "opt_records_to_remove.dat",
        include_iodine=False,
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
