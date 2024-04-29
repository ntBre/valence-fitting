#!/usr/bin/env bash
#SBATCH -J download
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

. ~/.bashrc

# Use the right conda environment
conda activate openff-nagl-test
conda env export > conda_env.yaml

# python download.py      \
#     --type      "optimization"                      \
#     --spec      "default"                           \
#     --output-directory  "input/optimization"        \
#     --dataset   "OpenFF Optimization Set 1"         \
#     --dataset   "SMIRNOFF Coverage Set 1"           \
#     --dataset   "OpenFF VEHICLe Set 1"              \
#     --dataset   "OpenFF Discrepancy Benchmark 1"    \
#     --dataset   "OpenFF Ehrman Informative Optimization v0.2"   \
#     --dataset   "Pfizer discrepancy optimization dataset 1"     \
#     --dataset   "FDA optimization dataset 1"                    \
#     --dataset   "Kinase Inhibitors: WBO Distributions"          \
#     --dataset   "OpenFF Gen 2 Opt Set 1 Roche"                  \
#     --dataset   "OpenFF Gen 2 Opt Set 2 Coverage"               \
#     --dataset   "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy"     \
#     --dataset   "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy" \
#     --dataset   "OpenFF Gen 2 Opt Set 5 Bayer"                  \
#     --dataset   "OpenFF Protein Fragments v1.0"                 \
#     --dataset   "OpenFF Protein Fragments v2.0"                 \
#     --dataset   "OpenFF Sandbox CHO PhAlkEthOH v1.0"            \
#     --dataset   "OpenFF Industry Benchmark Season 1 v1.1"       \
    # --dataset   "OpenFF Gen2 Optimization Dataset Protomers v1.0"  \
    # --dataset   "OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0" \
    # --dataset   "OpenFF Iodine Chemistry Optimization Dataset v1.0" \
    # --dataset   "XtalPi Shared Fragments OptimizationDataset v1.0"

# python download.py      \
#     --type      "torsiondrive"                      \
#     --spec      "default"                           \
#     --output-directory  "input/torsiondrive"        \
#     --dataset    "Fragment Stability Benchmark"       \
#     --dataset    "OpenFF Group1 Torsions"       \
#     --dataset    "SMIRNOFF Coverage Torsion Set 1"       \
#     --dataset    "OpenFF Substituted Phenyl Set 1"       \
#     --dataset    "Pfizer discrepancy torsion dataset 1"       \
#     --dataset    "TorsionDrive Paper"       \
#     --dataset    "OpenFF Primary Benchmark 1 Torsion Set"       \
#     --dataset    "OpenFF Primary Benchmark 2 Torsion Set"       \
#     --dataset    "OpenFF Group1 Torsions 2"       \
#     --dataset    "OpenFF Group1 Torsions 3"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 1 Roche"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 2 Coverage"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 5 Bayer"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 6 supplemental"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 1 Roche 2"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 2 Coverage 2"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 5 Bayer 2"       \
#     --dataset    "OpenFF Gen 2 Torsion Set 6 supplemental 2"       \
#     --dataset    "OpenFF Fragmenter Validation 1.0"       \
#     --dataset    "OpenFF DANCE 1 eMolecules t142 v1.0"       \
#     --dataset    "OpenFF Rowley Biaryl v1.0"       \
#     --dataset    "OpenFF-benchmark-ligand-fragments-v1.0"       \
#     --dataset    "OpenFF Protein Fragments TorsionDrives v1.0"       \
#     --dataset    "OpenFF WBO Conjugated Series v1.0"       \
#     --dataset    "OpenFF Amide Torsion Set v1.0"       \
#     --dataset    "OpenFF Gen3 Torsion Set v1.0"       \
#     --dataset    "OpenFF Aniline 2D Impropers v1.0"       \
#     --dataset    "OpenFF-benchmark-ligand-fragments-v2.0"       \
#     --dataset    "OpenFF Protein Dipeptide 2-D TorsionDrive v2.1"       \
#     --dataset    "OpenFF Protein Capped 1-mer Sidechains v1.3"       \
#     --dataset    "OpenFF Protein Capped 3-mer Backbones v1.0"       \
#     --dataset    "OpenFF multiplicity correction torsion drive data v1.1"       \
#     --dataset    "OpenFF Protein Capped 3-mer Omega v1.0"       \
#     --dataset    "XtalPi Shared Fragments TorsiondriveDataset v1.0"       \
#     --dataset    "OpenFF Torsion Coverage Supplement v1.0"       \
