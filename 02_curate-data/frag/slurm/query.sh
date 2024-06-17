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

date
hostname
echo jobid $SLURM_JOB_ID

source $HOME/.bashrc

conda activate fb-196-qcnew

ff=../../01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml

python query.py \
	   -f $ff \
	   -n $SLURM_CPUS_PER_TASK \
	   -t want.opt \
	   -x 'natoms:80' \
	   -x 'elements:Cl, P, Br, I, H, C, O, N, F, S' \
	   -x 'inchi:inchis.dat'

date
