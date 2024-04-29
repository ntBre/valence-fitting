import json
import os
import pathlib
import tempfile

import click
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from openff.toolkit import Molecule
from rdkit.Chem import Draw
from reportlab.graphics import renderPM
from svglib.svglib import svg2rlg

sns.set_context("talk")
mpl.rcParams["font.sans-serif"] = ["muli"]


def draw_grid_df(
    df,
    use_svg: bool = True,
    output_file: str = None,
    n_col: int = 4,
    n_page: int = 24,
    subImgSize=(300, 300),
):
    """
    Draw molecules

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the molecules to draw.
    use_svg : bool, optional
        Whether to use SVG format, by default True
    output_file : str, optional
        The output file to save the images, by default None.
        If None, the images are not saved.
        If there are more than `n_page` images,
        the images are saved in chunks.
    n_col : int, optional
        The number of columns in the grid, by default 4
    n_page : int, optional
        The number of images per page, by default 24
    subImgSize : tuple, optional
        The size of the subimages, by default (300, 300)
    """

    rdmols = []
    legends = []
    tags = []
    for _, row in df.iterrows():
        indices = row["atom_indices"]
        mol = Molecule.from_mapped_smiles(
            row["mapped_smiles"], allow_undefined_stereo=True
        )
        rdmol = mol.to_rdkit()
        for index in indices:
            atom = rdmol.GetAtomWithIdx(int(index))
            atom.SetProp("atomNote", str(index))
        rdmols.append(rdmol)
        indices_text = "-".join(list(map(str, indices)))
        legends.append(f"{row['torsiondrive_id']}: {indices_text}")
        tags.append(tuple(map(int, indices)))

    images = []
    for i in range(0, len(rdmols), n_page):
        j = i + n_page
        rdmols_chunk = rdmols[i:j]
        legends_chunk = legends[i:j]
        tags_chunk = tags[i:j]

        img = Draw.MolsToGridImage(
            rdmols_chunk,
            molsPerRow=n_col,
            legends=legends_chunk,
            subImgSize=subImgSize,
            # maxMols=n_page,
            highlightAtomLists=tags_chunk,
            returnPNG=False,
            useSVG=use_svg,
        )

        images.append(img)
    if output_file:
        output_file = pathlib.Path(output_file)
        output_file.parent.mkdir(exist_ok=True, parents=True)
        if not use_svg:
            images[0].save(
                output_file,
                append_images=images[1:],
                save_all=True,
                dpi=(300, 300),
            )
            print(f"Saved {output_file}")
        else:
            base_file, suffix = str(output_file).rsplit(".", maxsplit=1)
            for i, img in enumerate(images):
                file = f"{base_file}_{i}.{suffix}"

                with tempfile.TemporaryDirectory() as tempdir:
                    cwd = os.getcwd()
                    os.chdir(tempdir)

                    try:
                        data = img.data
                    except AttributeError:
                        data = img
                    with open("temp.svg", "w") as f:
                        f.write(data)
                    drawing = svg2rlg("temp.svg")
                os.chdir(cwd)
                renderPM.drawToFile(drawing, file, fmt="PNG")
                print(f"Saved {file}")

    return images


@click.command()
@click.option(
    "--forcefield",
    type=str,
    help="The forcefield to use.",
)
@click.option(
    "--parameter-id",
    type=str,
    help="The parameter id to plot.",
)
@click.option(
    "--output-directory",
    type=str,
    help="The directory to save the plots.",
    default="../images",
)
@click.option(
    "--qm-dataset",
    "qm_dataset_path",
    type=str,
    help="The path to the QM dataset.",
    default="datasets/qm/output/torsiondrive",
)
@click.option(
    "--mm-dataset",
    "mm_dataset_path",
    type=str,
    help="The path to the MM dataset.",
    default="datasets/mm/singlepoint-torsiondrive-datasets",
    required=False,
)
@click.option(
    "--parameter-ids-to-torsions",
    "parameter_ids_to_torsions_path",
    type=str,
    help="The path to the parameter id to torsion ids mapping.",
    default="parameter_id_to_torsion_ids.json",
)
@click.option(
    "--plot-molecules",
    is_flag=True,
    help="Whether to plot the molecules.",
)
@click.option(
    "--suffix",
    type=str,
    help="The suffix to append to the output file.",
    default="",
)
def main(
    forcefield: str,
    parameter_id: str,
    output_directory: str,
    qm_dataset_path: str = "datasets/qm/output/torsiondrive",
    mm_dataset_path: str = "datasets/mm/singlepoint-torsiondrive-datasets",
    parameter_ids_to_torsions_path: str = "parameter_id_to_torsion_ids.json",
    plot_molecules: bool = True,
    suffix: str = "",
):
    """
    Plot the QM and MM energies of torsion drives corresponding to
    a parameter ID. Optionally plot all molecules from the TorsionDrive
    as well.

    Dashed lines represent MM energies, solid lines represent QM energies.
    If dashed lines are not plotted, that means MM energies are not available.
    """

    import pyarrow.compute as pc
    import pyarrow.dataset as ds
    from openff.toolkit import ForceField

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

    if plot_molecules:
        draw_df = df.groupby("torsiondrive_id").first().reset_index()
        draw_grid_df(
            draw_df,
            output_file=output_directory / "molecules" / "molecules.png",
            subImgSize=(500, 500),
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
