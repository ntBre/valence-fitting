import click
from openff.qcsubmit.results import OptimizationResultCollection
from qcportal import PortalClient


@click.command()
@click.option("--output", "-o", default="datasets/core-opt.json")
@click.option(
    "--dataset",
    "-d",
    default="OpenFF multiplicity correction optimization set v1.0",
)
def main(output, dataset):
    td_datasets = [dataset]

    client = PortalClient("https://api.qcarchive.molssi.org:443/")
    dataset = OptimizationResultCollection.from_server(
        client=client,
        datasets=td_datasets,
        spec_name="default",
    )
    with open(output, "w") as out:
        out.write(dataset.json(indent=2))


if __name__ == "__main__":
    main()
