import json

import click
import numpy as np
import qcportal
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField

from curate_dataset import select_parameters


@click.group()
def cli():
    pass


@cli.command("select-td")
@click.option("--dataset", required=True)
@click.option("--forcefield", required=True)
@click.option("--output-smirks", required=True)
@click.option("--ring-torsions", required=True)
def select_td(dataset, forcefield, output_smirks, ring_torsions):
    dataset = TorsionDriveResultCollection.parse_file(dataset)

    ff = ForceField(
        forcefield,
        allow_cosmetic_attributes=True,
    )

    explicit_ring_torsions = np.loadtxt(ring_torsions, dtype=str)

    print("selecting td parameters")
    selected_parameters = select_parameters(
        dataset,
        ["ProperTorsions"],
        force_field=ff,
        explicit_ring_torsions=explicit_ring_torsions,
    )

    with open(output_smirks, "w") as f:
        json.dump(selected_parameters, f, indent=2)


@cli.command("select-impropers")
@click.option("--dataset", required=True)
@click.option("--forcefield", required=True)
@click.option("--output-smirks", required=True)
@click.option("--ring-torsions", required=True)
def select_impropers(dataset, forcefield, output_smirks, ring_torsions):
    dataset = TorsionDriveResultCollection.parse_file(dataset)

    ff = ForceField(
        forcefield,
        allow_cosmetic_attributes=True,
    )

    explicit_ring_torsions = np.loadtxt(ring_torsions, dtype=str)

    print("selecting improper parameters")
    selected_parameters = select_parameters(
        dataset,
        ["ImproperTorsions"],
        force_field=ff,
        explicit_ring_torsions=explicit_ring_torsions,
    )

    with open(output_smirks, "w") as f:
        json.dump(selected_parameters, f, indent=2)


@cli.command("select-opt")
@click.option("--dataset", required=True)
@click.option("--forcefield", required=True)
@click.option("--output-smirks", required=True)
def select_opt(dataset, forcefield, output_smirks):
    dataset = OptimizationResultCollection.parse_file(dataset)

    ff = ForceField(
        forcefield,
        allow_cosmetic_attributes=True,
    )

    print("selecting opt parameters")
    selected_parameters = select_parameters(
        dataset,
        ["Bonds", "Angles"],
        force_field=ff,
    )
    with open(output_smirks, "w") as file:
        json.dump(selected_parameters, file, indent=2)


if __name__ == "__main__":
    cli()
