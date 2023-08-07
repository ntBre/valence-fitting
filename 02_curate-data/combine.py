# Combine the sage 2.1 dataset with Pavan's

from openff.qcsubmit.results import (
    TorsionDriveResultCollection,
    OptimizationResultCollection,
)
import typing
import click

Output: typing.TypeAlias = typing.Union[
    "TorsionDriveResultCollection", "OptimizationResultCollection"
]


def combine_datasets(d1: Output, d2: Output) -> Output:
    "Combine `d1` and `d2` and filter the resulting collection for duplicates"
    key = list(d1.entries.keys())[0]

    all_entries = d1.entries[key] + d2.entries[key]

    output_type = type(d1)

    assert output_type == type(d2)

    # filter duplicates
    unique = {record.record_id: record for record in all_entries}
    new_dataset = (output_type)(entries={key: list(unique.values())})

    return new_dataset


@click.group()
def cli():
    pass


@cli.command("combine-td")
@click.option("--input-datasets", nargs=2, required=True)
@click.option("--output-dataset", required=True)
def combine_td(input_datasets, output_dataset):
    ds1, ds2 = input_datasets
    sage_td = TorsionDriveResultCollection.parse_file(ds1)
    pavan_td = TorsionDriveResultCollection.parse_file(ds2)

    print("combining td datasets")
    combined_td = combine_datasets(sage_td, pavan_td)

    with open(output_dataset, "w") as out:
        out.write(combined_td.json(indent=2))


@cli.command("combine-opt")
@click.option("--input-datasets", nargs=2, required=True)
@click.option("--output-dataset", required=True)
def combine_opt(input_datasets, output_dataset):
    ds1, ds2 = input_datasets
    sage_td = OptimizationResultCollection.parse_file(ds1)

    pavan_td = OptimizationResultCollection.parse_file(ds2)

    print("combining opt datasets")
    combined_td = combine_datasets(sage_td, pavan_td)

    with open(output_dataset, "w") as out:
        out.write(combined_td.json(indent=2))


if __name__ == "__main__":
    cli()
