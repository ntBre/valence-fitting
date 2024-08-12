import click
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import (
    ConformerRMSDFilter,
    ConnectivityFilter,
    ElementFilter,
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
    dataset = OptimizationResultCollection.parse_file(input)

    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

    dataset = dataset.filter(
        NoisyFilter(name="RecordStatusFilter"),
        RecordStatusFilter(status=RecordStatusEnum.complete),
        NoisyFilter(name="ConnectivityFilter"),
        ConnectivityFilter(tolerance=1.2),
        NoisyFilter(name="UnperceivableStereoFilter()"),
        UnperceivableStereoFilter(),
        NoisyFilter(name="ElementFilter"),
        ElementFilter(allowed_elements=elements),
        NoisyFilter(name="ConformerRMSDFilter"),
        ConformerRMSDFilter(max_conformers=12),  # default arg from filters
        NoisyFilter(name="ChargeCheckFilter"),
        ChargeCheckFilter(nprocs=nprocs, chunksize=chunksize),
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
