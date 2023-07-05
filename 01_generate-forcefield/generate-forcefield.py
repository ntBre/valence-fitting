import click
from dataclasses import dataclass
from openff.toolkit import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import (
    ParameterType,
)

TORSIONS = "ProperTorsions"


@dataclass
class Param:
    param: ParameterType
    after: str


def get_smirks(params):
    return [p.smirks for p in params]


def by_smirks(smirks):
    return dict(smirks=smirks)


@click.command()
@click.option(
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file.",
)
@click.option(
    "--force-field-name",
    type=str,
    default="openff-2.1.0.offxml",
    help="The name of the force field to download.",
    show_default=True,
)
def download_force_field(
    output_path: str,
    force_field_name: str = "openff-2.1.0.offxml",
):
    sage20 = ForceField("openff-2.0.0.offxml")
    tors20 = sage20.get_parameter_handler(TORSIONS)

    pavan = ForceField("force-field.offxml", allow_cosmetic_attributes=True)
    torspv = pavan.get_parameter_handler(TORSIONS)

    ids20 = get_smirks(tors20.parameters)
    idspv = get_smirks(torspv.parameters)

    set20 = set(ids20)
    setpv = set(idspv)

    # - only in pavan's set
    only_pavan = setpv.difference(set20)

    # gather the parameters only found in pavan's set
    new_params = []
    for i, param in enumerate(torspv.parameters):
        if param.smirks in only_pavan:
            new_params.append(Param(param, torspv.parameters[i - 1].smirks))

    # use 2.1 as the base for the output force field
    force_field = ForceField(force_field_name)
    h = force_field.get_parameter_handler(TORSIONS)

    # add the new parameters
    for p in new_params:
        # avoid duplicate parameter error
        if not h.get_parameter(by_smirks(p.param.smirks)):
            # ensure unique id
            if p.param.id in [p.id for p in h.parameters]:
                p.param.id += "x"
            h.add_parameter(parameter=p.param, after=p.after)

    force_field.to_file(output_path)


if __name__ == "__main__":
    download_force_field()
