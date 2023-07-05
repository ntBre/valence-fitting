#!/bin/bash


source ~/.bashrc
micromamba activate fb-195-tk-014-py310

python generate-forcefield.py                               \
    --force-field-name          "openff-2.1.0.offxml"       \
    --output                    "output/initial-force-field-openff-2.1.0.offxml"
