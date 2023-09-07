import click
from openff.qcsubmit.results import TorsionDriveResultCollection

from filters import filter_td_data


@click.command()
@click.option("--input")
@click.option("--output")
def main(input, output):
    dataset = TorsionDriveResultCollection.parse_file(input)
    dataset = filter_td_data(
        dataset,
        "td_records_to_remove.dat",
        include_iodine=False,
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
