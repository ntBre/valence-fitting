import click
import pprint
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


def by_id(_id):
    return dict(id=_id)


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
    initial_parameters = [p.id for p in h.parameters]

    # add the new parameters
    params_to_delete = []
    for p in new_params:
        # avoid duplicate parameter error
        if not h.get_parameter(by_smirks(p.param.smirks)):
            # ensure unique id
            if p.param.id in initial_parameters:
                params_to_delete.append(p.param.id)
                p.param.id += "x"
            h.add_parameter(parameter=p.param, after=p.after)

    # check if the right smirks were removed
    removed_in_multiplicity = set20 - setpv

    actually_removed = set(
        [h.get_parameter(by_id(p))[0].smirks for p in params_to_delete]
    )

    print("smirks removed in multiplicity but not actually removed: ")
    pprint.pprint(removed_in_multiplicity - actually_removed)

    print()

    print("smirks actually removed but shouldn't be:")
    pprint.pprint(actually_removed - removed_in_multiplicity)

    # remove the duplicate parameter ids from above
    lh = len(h.parameters)
    print("before=", lh)

    # pop from the end to avoid moving out of the part we're iterating over
    for i in reversed(range(lh)):
        p = h.parameters[i]
        if p.id in params_to_delete:
            h.parameters.pop(i)
    print("after=", len(h.parameters))

    force_field.to_file(output_path)


if __name__ == "__main__":
    download_force_field()
