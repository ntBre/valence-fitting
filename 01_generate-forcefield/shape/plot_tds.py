import json
import pathlib

import click
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
import seaborn as sns
from matplotlib import pyplot as plt
from openff.toolkit import ForceField

sns.set_context("talk")


@click.command()
@click.option("--forcefield")
@click.option("--parameter-id")
@click.option("--output-directory")
@click.option("--qm-dataset")
@click.option("--mm-dataset")
@click.option("--parameter-ids-to-torsions")
@click.option("--suffix")
def main(
    forcefield: str,
    parameter_id: str,
    output_directory: str,
    qm_dataset_path: str = "datasets/qm/output/torsiondrive",
    mm_dataset_path: str = "datasets/mm/singlepoint-torsiondrive-datasets",
    parameter_ids_to_torsions_path: str = "parameter_id_to_torsion_ids.json",
    suffix: str = "",
):
    """
    Plot the QM and MM energies of torsion drives corresponding to
    a parameter ID. Optionally plot all molecules from the TorsionDrive
    as well.

    Dashed lines represent MM energies, solid lines represent QM energies.
    If dashed lines are not plotted, that means MM energies are not available.
    """

    qm_dataset = ds.dataset(qm_dataset_path)
    mm_dataset = ds.dataset(mm_dataset_path)
    if "forcefield" in mm_dataset.schema.names:
        mm_dataset = mm_dataset.filter(pc.field("forcefield") == forcefield)

    with open(parameter_ids_to_torsions_path, "r") as f:
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
        "mapped_smiles",
        "dihedral",
        "qcarchive_id",
    ]
    df = subset.to_table(columns=cols).to_pandas()
    df["atom_indices"] = [tuple(x) for x in df["dihedral"]]

    output_directory = (
        pathlib.Path(output_directory)
        / pathlib.Path(forcefield).stem
        / parameter_id
    )

    # get MM energies
    mm_expression = pc.field("qcarchive_id").isin(df.qcarchive_id.values)
    mm_subset = mm_dataset.filter(mm_expression)
    mm_energies = mm_subset.to_table(columns=["qcarchive_id", "mm_energy"])
    mm_df = mm_energies.to_pandas().sort_values("qcarchive_id")

    if len(mm_df) == 0:
        print(f"{parameter_id} has no minimized geometries")
        df["mm_energy"] = np.nan
    else:
        df = df.sort_values("qcarchive_id")
        df = df.merge(
            mm_df,
            left_on=["qcarchive_id"],
            right_on=["qcarchive_id"],
            how="left",
        )

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
        relative_qm = subdf["energy"] - subdf["energy"].values[lowest_index]
        subdf["relative_qm_energy"] = relative_qm
        relative_mm = (
            subdf["mm_energy"] - subdf["mm_energy"].values[lowest_index]
        )
        subdf["relative_mm_energy"] = relative_mm
        subdfs.append(subdf)
    df = pd.concat(subdfs)
    df = df.melt(
        id_vars=[
            x
            for x in df.columns
            if x not in ["relative_qm_energy", "relative_mm_energy"]
        ],
        value_vars=["relative_qm_energy", "relative_mm_energy"],
        var_name="Type",
        value_name="relative_energy",
    )
    df = df.sort_values(by=["torsiondrive_id", "grid_id"])
    df["TorsionDrive ID"] = [str(x) for x in df.torsiondrive_id.values]

    plt.clf()
    g = sns.FacetGrid(data=df, aspect=1.4, height=4, hue="TorsionDrive ID")
    g.map_dataframe(sns.lineplot, "grid_id", "relative_energy", style="Type")
    title = f"{parameter_id}\n" + smirks.encode("unicode_escape").decode(
        "utf-8"
    ).replace("$", "\\$")
    ax = list(g.axes.flatten())[0]
    ax.set_title(title)
    ax.set_xlabel("Angle (°)")
    ax.set_ylabel("Relative energy\n[kcal/mol]")
    plt.tight_layout()

    g.add_legend()
    filename = output_directory / f"mm-torsion-energies{suffix}.png"
    g.savefig(filename, dpi=300)
    print(f"Saved to {filename}")


if __name__ == "__main__":
    main()
