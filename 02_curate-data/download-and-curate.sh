#!/bin/bash

set -e

init_ff="../01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml"

python curate_dataset.py download-opt                                                  \
    --core-opt-dataset          "OpenFF multiplicity correction optimization set v1.0" \
    --initial-forcefield        $init_ff                                               \
    --max-opt-conformers        12                                                     \
    --n-processes               8                                                      \
    --min-record-coverage       1                                                      \
    --output                    "output/pavan-opt-training-set.json"                   \
    --output-parameter-smirks   "output/pavan-opt-smirks.json"                         \
    --verbose

python curate_dataset.py download-td                                                     \
    --core-td-dataset           "OpenFF multiplicity correction torsion drive data v1.1" \
    --initial-forcefield        $init_ff                                                 \
    --n-processes               8                                                        \
    --min-record-coverage       1                                                        \
    --output                    "output/pavan-td-training-set.json"                      \
    --output-parameter-smirks   "output/pavan-td-torsion-smirks.json"                    \
    --verbose
