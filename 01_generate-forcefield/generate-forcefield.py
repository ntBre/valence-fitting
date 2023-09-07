import copy
import re
from dataclasses import dataclass

import click
from openff.toolkit import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import ParameterType

TORSIONS = "ProperTorsions"


@dataclass
class Param:
    param: ParameterType
    after: str


def get_smirks(params):
    return [p.smirks for p in params]


def get_ids(params):
    return [p.id for p in params]


def by_smirks(smirks):
    return dict(smirks=smirks)


def by_id(_id):
    return dict(id=_id)


def numeric(param):
    "Sort key function for parameter ids, numerically then by suffix"
    m = re.match(r"t(\d+)([a-z]+)?", param.id)
    num, let = m.groups()
    if not let:
        let = "A"
    return int(num), let


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

    sids = set(get_ids(tors20))
    pids = set(get_ids(torspv))

    # these are the torsion parameter ids removed in the torsion multiplicity
    # work, presumably to be replaced by more specific variants
    removed_by_pavan: set[str] = sids - pids

    # Sage 2.1.0 force field to work on
    ret = ForceField(force_field_name)
    h = ret.get_parameter_handler(TORSIONS)

    print(f"initial torsions: {len(h.parameters)}")

    # parameter handlers are so annoying to work with. convert these ids to
    # indices, sort them in descending order, and pop them off the parameter
    # list
    indices = [
        i for i, p in enumerate(h.parameters) if p.id in removed_by_pavan
    ]
    for i in sorted(indices, reverse=True):
        h.parameters.pop(i)

    print(f"removed by pavan: {len(h.parameters)}")

    # at this point we've deleted the parameters pavan deleted from sage 2.0
    # and can begin adding the parameters he added
    added_by_pavan: set[str] = pids - sids
    for pid in added_by_pavan:
        param = torspv.get_parameter(by_id(pid))[0]
        # add an x to our parameters so they sort after the existing one
        if h.get_parameter(by_id(pid)):
            param.id += "x"
        h.add_parameter(parameter=param)

    print(f"added by pavan: {len(h.parameters)}")

    params = copy.deepcopy(h.parameters)
    h.parameters.clear()

    params = sorted(params, key=numeric)
    for p in params:
        h.add_parameter(parameter=p)

    print(f"final: {len(h.parameters)}")

    ret.to_file(output_path)


if __name__ == "__main__":
    download_force_field()
