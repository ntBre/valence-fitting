import click
from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    ElementFilter,
    HydrogenBondFilter,
    RecordStatusEnum,
    RecordStatusFilter,
    UnperceivableStereoFilter,
)

from filters import ChargeCheckFilter


@click.command()
@click.option("--input")
@click.option("--output")
def main(input, output):
    "Filter a TorsionDrive dataset"

    dataset = TorsionDriveResultCollection.parse_file(input)

    records_to_remove = []

    key = list(dataset.entries.keys())[0]

    # filter out entries to remove
    dataset.entries[key] = [
        entry
        for entry in dataset.entries[key]
        if entry.record_id not in records_to_remove
    ]

    # in a number of datasets the iodine-containing molecules
    # were tainted due to an auxiliary basis set issue
    # This has since been resolved and entries have been recomputed
    # in new datasets, but we still need to filter the old ones
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]

    # filter out other unsuitable entries
    dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        HydrogenBondFilter(method="baker-hubbard"),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        ElementFilter(allowed_elements=elements),
        ChargeCheckFilter(),
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
