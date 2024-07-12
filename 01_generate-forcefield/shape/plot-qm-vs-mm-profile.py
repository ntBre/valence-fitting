import json
import os
import pathlib
import tempfile

import click
import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds
import seaborn as sns
import tqdm
from matplotlib import pyplot as plt

sns.set_context("talk")

OPENFF_BLUE = "#015480"
OPENFF_LIGHT_BLUE = "#2F9ED2"
OPENFF_ORANGE = "#F08521"
OPENFF_RED = "#F03A21"
OPENFF_GRAY = "#3E424A"

COLORS = {
    "blue": OPENFF_BLUE,
    "cyan": OPENFF_LIGHT_BLUE,
    "orange": OPENFF_ORANGE,
    "red": OPENFF_RED,
    "gray": OPENFF_GRAY,
}


def draw_single(df, torsiondrive_id: int, width=300, height=300):
    from cairosvg import svg2png
    from matplotlib import pyplot as plt
    from openff.toolkit import Molecule
    from rdkit.Chem.Draw import rdMolDraw2D

    subdf = df[df.torsiondrive_id == torsiondrive_id]
    mol = Molecule.from_mapped_smiles(
        subdf["mapped_smiles"].values[0], allow_undefined_stereo=True
    )
    rdmol = mol.to_rdkit()
    dihedral = list(map(int, subdf["dihedral"].values[0]))
    for index in dihedral:
        atom = rdmol.GetAtomWithIdx(int(index))
        atom.SetProp("atomNote", str(index))
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    # options.baseFontSize = 1
    drawer.DrawMolecule(rdmol, highlightAtoms=dihedral)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    with tempfile.TemporaryDirectory() as tempdir:
        cwd = os.getcwd()
        os.chdir(tempdir)
        svg2png(bytestring=svg, write_to="tmp.png", scale=10)
        png = plt.imread("tmp.png")
        os.chdir(cwd)

    return png


def plot_mm_vs_qm_profile(
    torsiondrive_id: int,
    mm_dataset,
    qm_dataset,
    output_directory=None,
    with_rmsds: bool = True,
    suffix: str = "",
):
    subset = qm_dataset.filter(pc.field("torsiondrive_id") == torsiondrive_id)
    df = subset.to_table().to_pandas()
    columns = ["qcarchive_id", "qm_energy", "mm_energy"]
    if with_rmsds:
        columns.append("RMSD_AA")
    geometry_df = (
        mm_dataset.filter(
            pc.field("qcarchive_id").isin(df.qcarchive_id.values)
        )
        .to_table(columns=columns)
        .to_pandas()
    )

    joined = df.merge(
        geometry_df,
        left_on=["qcarchive_id"],
        right_on=["qcarchive_id"],
        how="inner",
    )
    if not len(joined):
        print(f"{torsiondrive_id} has no minimized geometries: {len(df)} QM")
        return
    min_index = np.argmin(joined.qm_energy.values)
    joined["QM"] = joined.qm_energy - joined.qm_energy.values[min_index]
    joined["MM"] = joined.mm_energy - joined.mm_energy.values[min_index]
    joined = joined.sort_values("grid_id")

    fig, (ax1, imgax) = plt.subplots(figsize=(8, 4), ncols=2)
    color1 = "#015480"
    color1b = "#2F9ED2"
    ax1.set_xlabel("Angle")
    ax1.set_ylabel("Energy [kcal/mol]", color=color1)
    ax1.plot(joined.grid_id, joined.QM, color=color1, label="QM")
    ax1.plot(joined.grid_id, joined.MM, color=color1b, label="MM")
    ax1.tick_params(axis="y", labelcolor=color1)
    ax1.legend()

    if with_rmsds:
        ax2 = ax1.twinx()
        color2 = "#F03A21"
        ax2.set_ylabel("RMSD [Å]", color=color2)
        ax2.plot(joined.grid_id, joined.RMSD_AA, color=color2, label="RMSD")
        ax2.axhline(0.4, color=color2, ls="--", lw=1)
        ax2.tick_params(axis="y", labelcolor=color2)

        ax2.set_title(f"{torsiondrive_id}")
    else:
        ax1.set_title(f"{torsiondrive_id}")

    # png = draw_single(df, torsiondrive_id)
    # imgax.imshow(png, rasterized=True)
    imgax.set_xticks([])
    imgax.set_yticks([])
    imgax.spines["left"].set_visible(False)
    imgax.spines["right"].set_visible(False)
    imgax.spines["top"].set_visible(False)
    imgax.spines["bottom"].set_visible(False)

    fig.tight_layout()

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    basename = "qm-vs-mm"
    if with_rmsds:
        basename += "-rmsd"
    imgfile = output_directory / f"{basename}{suffix}.png"
    plt.savefig(imgfile, dpi=300)
    print(f"Saved to {imgfile}")


@click.command()
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
)
@click.option(
    "--forcefield",
    type=str,
    help="The forcefield to use.",
    default="tm-2.2.offxml",
)
@click.option(
    "--with-rmsds",
    is_flag=True,
    help="Whether to plot RMSDs.",
)
@click.option(
    "--parameter-ids-to-torsions",
    "parameter_ids_to_torsions_path",
    type=str,
    help="The path to the parameter id to torsion ids mapping.",
    default="parameter_id_to_torsion_ids.json",
)
@click.option(
    "--suffix",
    type=str,
    help="The suffix to append to the output file.",
    default="",
)
def plot_all(
    parameter_id: str,
    with_rmsds: bool = False,
    output_directory: str = "../images",
    qm_dataset_path: str = "datasets/qm/output/torsiondrive",
    mm_dataset_path: str = "datasets/mm/singlepoint-torsiondrive-datasets",
    forcefield: str = "tm-2.2.offxml",
    parameter_ids_to_torsions_path: str = "parameter_id_to_torsion_ids.json",
    suffix: str = "",
):
    """
    Plot QM vs MM profile with optional RMSDs if using minimized geometries.

    Outputs are saved in the following structure:
    output_directory/forcefield/parameter_id/qm-vs-mm/torsion_id.png
    """
    qm_dataset = ds.dataset(qm_dataset_path)
    mm_dataset = ds.dataset(mm_dataset_path)
    if "forcefield" in mm_dataset.schema.names:
        mm_dataset = mm_dataset.filter(pc.field("forcefield") == forcefield)

    with open(parameter_ids_to_torsions_path, "r") as f:
        parameter_id_to_torsion_ids = json.load(f)

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)
    ff_name = pathlib.Path(forcefield).stem

    torsion_ids = parameter_id_to_torsion_ids[parameter_id]
    for torsion_id in tqdm.tqdm(torsion_ids):
        plot_mm_vs_qm_profile(
            torsion_id,
            mm_dataset,
            qm_dataset,
            output_directory=(
                output_directory / ff_name / parameter_id / str(torsion_id)
            ),
            with_rmsds=with_rmsds,
            suffix=suffix,
        )


if __name__ == "__main__":
    plot_all()
