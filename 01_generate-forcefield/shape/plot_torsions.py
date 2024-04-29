import pathlib

import click
import numpy as np
import tqdm
from matplotlib import pyplot as plt
from openff.toolkit import ForceField
from openff.units import unit


def calc_torsion_energy(angle, parameter):
    angle = (angle * unit.degrees).m_as(unit.radians)
    total = 0 * unit.kilojoules_per_mole
    for k, phase, periodicity in zip(
        parameter.k, parameter.phase, parameter.periodicity
    ):
        phase = phase.m_as(unit.radians)
        total += k * (1 + np.cos(periodicity * angle - phase))
    return total


@click.command()
@click.option("--forcefield", "-f")
@click.option("--outdir", "-o")
def plot_ff_torsions(forcefield, outdir):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(exist_ok=True)
    ff = ForceField(forcefield, allow_cosmetic_attributes=True)
    handler = ff.get_parameter_handler("ProperTorsions")
    for parameter in tqdm.tqdm(handler.parameters):
        fig, ax = plt.subplots(figsize=(4, 3))
        xs = np.linspace(-180, 180, 360)
        ys = [
            calc_torsion_energy(x, parameter).m_as(unit.kilocalories_per_mole)
            for x in xs
        ]
        ax.plot(xs, ys)
        ax.set_title(parameter.id)
        ax.set_ylabel("Energy\n[kcal/mol]")
        plt.tight_layout()
        filename = outdir / f"{parameter.id}.png"

        plt.savefig(filename, dpi=300)
        plt.close()


if __name__ == "__main__":
    plot_ff_torsions()
