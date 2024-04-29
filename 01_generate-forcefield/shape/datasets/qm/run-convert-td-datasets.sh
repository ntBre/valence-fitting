#!/usr/bin/env bash
#SBATCH -J convert-datasets
#SBATCH --array=1-38
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A-%a.out

. ~/.bashrc

# Use the right conda environment
conda activate openff-nagl-test
conda env export > conda_env.yaml


TD_DATASET=(   \
    "Fragment-Stability-Benchmark"			\
    "OpenFF-Group1-Torsions-2"			\
    "OpenFF-Amide-Torsion-Set-v1.0"			\
    "OpenFF-Group1-Torsions-3"			\
    "OpenFF-Aniline-2D-Impropers-v1.0"			\
    "OpenFF-Group1-Torsions"			\
    "OpenFF-benchmark-ligand-fragments-v1.0"			\
    "OpenFF-multiplicity-correction-torsion-drive-data-v1.1"			\
    "OpenFF-benchmark-ligand-fragments-v2.0"			\
    "OpenFF-Primary-Benchmark-1-Torsion-Set"			\
    "OpenFF-DANCE-1-eMolecules-t142-v1.0"			\
    "OpenFF-Primary-Benchmark-2-Torsion-Set"			\
    "OpenFF-Fragmenter-Validation-1.0"			\
    "OpenFF-Protein-Capped-1-mer-Sidechains-v1.3"			\
    "OpenFF-Gen-2-Torsion-Set-1-Roche-2"			\
    "OpenFF-Protein-Capped-3-mer-Backbones-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-1-Roche"			\
    "OpenFF-Protein-Capped-3-mer-Omega-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-2-Coverage-2"			\
    "OpenFF-Protein-Dipeptide-2-D-TorsionDrive-v2.1"			\
    "OpenFF-Gen-2-Torsion-Set-2-Coverage"			\
    "OpenFF-Protein-Fragments-TorsionDrives-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy-2"			\
    "OpenFF-Rowley-Biaryl-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-3-Pfizer-Discrepancy"			\
    "OpenFF-Substituted-Phenyl-Set-1"			\
    "OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy-2"			\
    "OpenFF-Torsion-Coverage-Supplement-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-4-eMolecules-Discrepancy"			\
    "OpenFF-WBO-Conjugated-Series-v1.0"			\
    "OpenFF-Gen-2-Torsion-Set-5-Bayer-2"			\
    "Pfizer-discrepancy-torsion-dataset-1"			\
    "OpenFF-Gen-2-Torsion-Set-5-Bayer"			\
    "SMIRNOFF-Coverage-Torsion-Set-1"			\
    "OpenFF-Gen-2-Torsion-Set-6-supplemental-2"			\
    "TorsionDrive-Paper"			\
    "OpenFF-Gen-2-Torsion-Set-6-supplemental"			\
    "XtalPi-Shared-Fragments-TorsiondriveDataset-v1.0"			\
    "OpenFF-Gen3-Torsion-Set-v1.0"
)

TD_DATASET_NAME="${TD_DATASET[$SLURM_ARRAY_TASK_ID]}"

echo $TD_DATASET_NAME

python convert-td-to-dataset.py            \
    --input     "input/torsiondrive/${TD_DATASET_NAME}.json"   \
    --name      $TD_DATASET_NAME                               \
    --output    "output/torsiondrive/${TD_DATASET_NAME}"
