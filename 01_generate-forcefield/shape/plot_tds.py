import json
import pathlib

import click
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
import seaborn as sns
from matplotlib import pyplot as plt
from openff.toolkit import ForceField
from openff.units import unit

from plot_ff_torsions import calc_torsion_energy

sns.set_context("talk")


@click.command()
@click.option("--forcefield", "-f")
@click.option("--parameter-id", "-p")
@click.option("--output-directory", default="output")
@click.option("--qm-dataset", "-q")
@click.option("--parameter-ids-to-torsions", "-l")
@click.option("--suffix", default="")
def main(
    forcefield,
    parameter_id,
    output_directory,
    qm_dataset,
    parameter_ids_to_torsions,
    suffix,
):
    qm_dataset = ds.dataset(qm_dataset)

    with open(parameter_ids_to_torsions, "r") as f:
        parameter_id_to_torsion_ids = json.load(f)

    torsiondrive_ids = parameter_id_to_torsion_ids[parameter_id]
    if len(torsiondrive_ids) == 0:
        print(f"{parameter_id} has no torsiondrives")
        return

    ff = ForceField(forcefield, allow_cosmetic_attributes=True)
    handler = ff.get_parameter_handler("ProperTorsions")
    parameter = handler.get_parameter({"id": parameter_id})[0]
    smirks = parameter.smirks

    expression = pc.field("torsiondrive_id").isin(torsiondrive_ids)
    subset = qm_dataset.filter(expression)
    cols = [
        "torsiondrive_id",
        "grid_id",
        "energy",
        "dihedral",
        "qcarchive_id",
    ]
    df = subset.to_table(columns=cols).to_pandas()
    df["atom_indices"] = [tuple(x) for x in df["dihedral"]]
    df["ff_value"] = [
        calc_torsion_energy(x, parameter).m_as(unit.kilocalories_per_mole)
        for x in df.grid_id
    ]
    df["ff_value"] -= df["ff_value"].min()

    # remove accidental duplicates
    df = (
        df.groupby(
            by=["torsiondrive_id", "grid_id", "atom_indices", "qcarchive_id"]
        )
        .first()
        .reset_index()
    )

    # calculate relative energy
    subdfs = []
    for _, subdf in df.groupby("torsiondrive_id"):
        lowest_index = subdf["energy"].idxmin()
        relative_qm = subdf["energy"] - subdf["energy"][lowest_index]
        subdf["relative_qm_energy"] = relative_qm
        subdfs.append(subdf)
    df = pd.concat(subdfs)
    df = df.melt(
        id_vars=[x for x in df.columns if x not in ["relative_qm_energy"]],
        value_vars=["relative_qm_energy"],
        var_name="Type",
        value_name="relative_energy",
    )
    df = df.sort_values(by=["torsiondrive_id", "grid_id"])
    mn, mx = df.relative_energy.min(), df.relative_energy.max()
    df["ff_value"] *= mx - mn

    g = sns.FacetGrid(data=df, aspect=1.4, height=4, hue="torsiondrive_id")
    g.map_dataframe(sns.lineplot, "grid_id", "relative_energy", style="Type")
    g.map_dataframe(sns.lineplot, "grid_id", "ff_value", linestyle="--")

    ax = list(g.axes.flatten())[0]
    ax.set_title(f"{parameter_id}\n{smirks}")
    ax.set_xlabel("Angle (°)")
    ax.set_ylabel("Relative energy\n(kcal/mol)")
    plt.tight_layout()

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True)
    filename = output_directory / f"{parameter_id}.png"
    g.savefig(filename, dpi=300)
    print(f"Saved to {filename}")


if __name__ == "__main__":
    main()
