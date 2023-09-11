import click
from openff.qcsubmit.results import TorsionDriveResultCollection

from filters import filter_td_data


@click.command()
@click.option("--input")
@click.option("--output")
@click.option("--records-to-remove", default="td_records_to_remove.dat")
def main(input, output, records_to_remove):
    dataset = TorsionDriveResultCollection.parse_file(input)
    dataset = filter_td_data(
        dataset,
        records_to_remove,
        include_iodine=False,
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
