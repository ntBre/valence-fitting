#!/bin/bash
#SBATCH -J query
#SBATCH -p standard
#SBATCH -t 144:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=140gb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=bwestbr1@uci.edu
#SBATCH --constraint=fastscratch
#SBATCH -o logs/query.out

source $HOME/.bashrc

conda activate fb-196-qcnew

python query.py \
	   -n $SLURM_CPUS_PER_TASK \
	   -x 'natoms:80' \
	   -x 'elements:Cl, P, Br, I, H, C, O, N, F, S'
