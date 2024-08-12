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

from filters import ChargeCheckFilter, NoisyFilter


@click.command()
@click.option("--input")
@click.option("--output")
@click.option("--nprocs", "-n", default=1)
@click.option("--chunksize", "-c", default=1)
def main(input, output, nprocs, chunksize):
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

    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

    # filter out other unsuitable entries
    dataset = dataset.filter(
        NoisyFilter(name="RecordStatusFilter"),
        RecordStatusFilter(status=RecordStatusEnum.complete),
        NoisyFilter(name="HydrogenBondFilter"),
        HydrogenBondFilter(method="baker-hubbard"),
        NoisyFilter(name="ConnectivityFilter"),
        ConnectivityFilter(tolerance=1.2),
        NoisyFilter(name="UnperceivableStereoFilter()"),
        UnperceivableStereoFilter(),
        NoisyFilter(name="ElementFilter"),
        ElementFilter(allowed_elements=elements),
        NoisyFilter(name="ChargeCheckFilter"),
        ChargeCheckFilter(nprocs=nprocs, chunksize=chunksize),
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
