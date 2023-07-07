#!/bin/bash

python create-fb-inputs.py                                                                          \
    --tag                       "fb-fit"                                                            \
    --optimization-dataset      "../02_curate-data/output/combined-opt.json"           \
    --torsion-dataset           "../02_curate-data/output/combined-td.json"                \
    --forcefield                "../03_generate-initial-ff/output/initial-force-field-msm.offxml"   \
    --valence-to-optimize       "../02_curate-data/output/combined-opt-smirks.json"             \
    --torsions-to-optimize      "../02_curate-data/output/combined-td-smirks.json"             \
    --smiles-to-exclude         "smiles-to-exclude.dat"                                             \
    --smarts-to-exclude         "smarts-to-exclude.dat"                                             \
    --max-iterations            1                                                                   \
    --port                      55387                                                               \
    --output-directory          "output"                                                            \
    --verbose

