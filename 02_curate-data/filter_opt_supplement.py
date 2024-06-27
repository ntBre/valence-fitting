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

from filters import ChargeCheckFilter


@click.command()
@click.option("--input")
@click.option("--output")
def main(input, output):
    dataset = OptimizationResultCollection.parse_file(input)

    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

    dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        ElementFilter(allowed_elements=elements),
        ConformerRMSDFilter(max_conformers=12),  # default arg from filters
        ChargeCheckFilter(),
    )

    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
