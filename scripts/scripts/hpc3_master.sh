#!/bin/bash
#SBATCH -J valence-fit
#SBATCH -p standard
#SBATCH -t 144:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000mb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=bwestbr1@uci.edu
#SBATCH --constraint=fastscratch
#SBATCH -o master.out
#SBATCH -e master.err

conda_env=fb-196

compressed_env=/dfs4/dmobley-lab/$USER/envs/$conda_env.tar.gz

TMPDIR=/tmp/$USER/$SLURM_JOB_ID

rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR

source $HOME/.bashrc

cp $compressed_env .
mkdir -p $conda_env
tar xzf $compressed_env -C $conda_env
mamba activate $conda_env

echo $CONDA_PREFIX > $SLURM_SUBMIT_DIR/env.path

scp -C  $SLURM_SUBMIT_DIR/optimize.in     $TMPDIR
scp -C  $SLURM_SUBMIT_DIR/targets.tar.gz  $TMPDIR
scp -Cr $SLURM_SUBMIT_DIR/forcefield      $TMPDIR

tar -xzf targets.tar.gz

datadir=$(pwd)
echo $(hostname) > $SLURM_SUBMIT_DIR/host

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

if ForceBalance.py optimize.in ; then
   tar -czvf optimize.tmp.tar.gz optimize.tmp
   rsync  -avzIi --exclude="optimize.tmp" --exclude="optimize.bak" \
	  --exclude="fb_193*" \
	  --exclude="targets*" $TMPDIR/* $SLURM_SUBMIT_DIR > copy.log
   rm -rf $TMPDIR
fi

echo "All done"
