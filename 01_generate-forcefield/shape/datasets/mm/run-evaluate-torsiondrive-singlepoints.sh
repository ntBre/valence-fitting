#!/usr/bin/env bash
#SBATCH -J singlepoint-torsions
#SBATCH -p free
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate ib-dev

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"


# VERSION="2.2.0-rc1"
VERSION="tm-2.2"

echo $VERSION

python evaluate-torsiondrive-singlepoints.py                          \
        --n-workers                     300                     \
        --worker-type                   "slurm"                 \
        --batch-size                    500                     \
        --memory                        4                       \
        --walltime                      480                     \
        --queue                         "free"                  \
        --conda-environment             "ib-dev"                \
    --input         "/dfs9/dmobley-lab/lilyw7/datasets/qm/output/torsiondrive" \
    --force-field   "${VERSION}.offxml"          \
    --output        "singlepoint-torsiondrive-datasets/${VERSION}"     \

