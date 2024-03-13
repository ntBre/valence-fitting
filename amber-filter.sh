#!/bin/bash
#SBATCH -J amber-filter
#SBATCH -p standard
#SBATCH -t 144:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=bwestbr1@uci.edu
#SBATCH --constraint=fastscratch
#SBATCH -o af.out

source $HOME/.bashrc
export OE_LICENSE=""

conda activate fb-196

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

make 02_curate-data/sage/amber-filtered-td.json
make 02_curate-data/sage/amber-filtered-opt.json
