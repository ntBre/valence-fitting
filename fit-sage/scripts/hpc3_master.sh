#!/bin/bash
#SBATCH -J valence-fit
#SBATCH -p standard
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000mb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=MY_USERNAME@uci.edu
#SBATCH --constraint=fastscratch
#SBATCH -o master.out
#SBATCH -e master.err

rm -rf /tmp/$SLURM_JOB_NAME
source $HOME/.bashrc

mamba activate valence-fitting

rsync  -avzIi  $SLURM_SUBMIT_DIR/optimize.in  $SLURM_TMPDIR/$SLURM_JOB_NAME
rsync  -avzIi  $SLURM_SUBMIT_DIR/targets.tar.gz  $SLURM_TMPDIR/$SLURM_JOB_NAME
rsync  -avzIi  $SLURM_SUBMIT_DIR/forcefield  $SLURM_TMPDIR/$SLURM_JOB_NAME

tar -xzf targets.tar.gz

datadir=$(pwd)
mkdir -p $SLURM_SUBMIT_DIR/worker_logs
echo $(hostname) > $SLURM_SUBMIT_DIR/host

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

if ForceBalance.py optimize.in ; then
   tar -czvf optimize.tmp.tar.gz optimize.tmp
   rsync  -avzIi --exclude="optimize.tmp" --exclude="optimize.bak" \
	  --exclude="fb_193*" --exclude="targets*" $TMPDIR/* \
	  $SLURM_SUBMIT_DIR > copy.log
   rm -rf $TMPDIR
fi

echo "All done"
