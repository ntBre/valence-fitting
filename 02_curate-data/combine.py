# Combine the sage 2.1 dataset with Pavan's

import typing

import click
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)

Output: typing.TypeAlias = typing.Union[
    "TorsionDriveResultCollection", "OptimizationResultCollection"
]


def combine_datasets(datasets) -> Output:
    "Combine `d1` and `d2` and filter the resulting collection for duplicates"
    assert len(datasets) > 0
    assert len(datasets[0].entries.keys()) == 1

    key = list(datasets[0].entries.keys())[0]

    output_type = type(datasets[0])

    all_entries = []
    for ds in datasets:
        assert output_type == type(ds)
        all_entries += ds.entries[key]

    # filter duplicates
    unique = {record.record_id: record for record in all_entries}
    new_dataset = (output_type)(entries={key: list(unique.values())})

    return new_dataset


@click.group()
def cli():
    pass


def common(cls, input_datasets, output_dataset):
    datasets = [cls.parse_file(ds) for ds in input_datasets]
    combined_td = combine_datasets(datasets)
    with open(output_dataset, "w") as out:
        out.write(combined_td.json(indent=2))


@cli.command("combine-td")
@click.option("--input-datasets", multiple=True, required=True)
@click.option("--output-dataset", required=True)
def combine_td(input_datasets, output_dataset):
    common(TorsionDriveResultCollection, input_datasets, output_dataset)


@cli.command("combine-opt")
@click.option("--input-datasets", multiple=True, required=True)
@click.option("--output-dataset", required=True)
def combine_opt(input_datasets, output_dataset):
    common(OptimizationResultCollection, input_datasets, output_dataset)


if __name__ == "__main__":
    cli()
