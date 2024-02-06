import click
from openff.qcsubmit.results.filters import ResultRecordFilter
from openff.toolkit.utils.exceptions import (
    ChargeCalculationError,
    ConformerGenerationError,
)
from openff.toolkit.utils.toolkits import AmberToolsToolkitWrapper
from vflib import Molecules, load_dataset

counter = 0


class ChargeCheckFilter(ResultRecordFilter):
    def _filter_function(self, result, record, molecule) -> bool:
        global counter
        if counter % 100 == 0:
            print(counter)
        counter += 1
        try:
            AmberToolsToolkitWrapper().assign_partial_charges(
                molecule, partial_charge_method="am1bcc"
            )
            # from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
            # from openff.units import unit
            # molecule.generate_conformers(
            #     n_conformers=1,
            #     rms_cutoff=0.25 * unit.angstrom,
            #     toolkit_registry=RDKitToolkitWrapper(),
            # )
        except (ChargeCalculationError, ConformerGenerationError):
            return False
        else:
            return True


@click.command()
@click.option("--input")
@click.option("--output")
def main(input, output):
    dataset = load_dataset(input)
    init = len(Molecules(dataset))
    print(f"loaded {init} entries from {input}")
    dataset = dataset.filter(ChargeCheckFilter())
    with open(output, "w") as out:
        out.write(dataset.json(indent=2))
    final = len(Molecules(dataset))
    print(f"filtered {init} entries in {input} to {final}")


if __name__ == "__main__":
    main()
