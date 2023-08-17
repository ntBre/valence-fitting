#!/bin/bash

host=$(sed 1q host)
port=$(awk '/port/ {print $NF}' optimize.in)

USERNAME=$(whoami)
export SLURM_TMPDIR=/tmp
export MYTMPDIR="/tmp/${USERNAME}"
export TMPDIR=$SLURM_TMPDIR/$SLURM_JOB_NAME

worker_num=$(squeue -u $USER | grep wq -c)
ncpus=10

echo host = $host
echo port = $port
echo $ncpus cpus requested
echo submitting worker $worker_num

#export OE_LICENSE=

cmd=$(mktemp)
cat << EOF > $cmd
#!/usr/bin/env bash
#SBATCH -J wq-$port
#SBATCH -p free
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=${ncpus}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
# SBATCH --array=1-100
#SBATCH --account dmobley_lab
# SBATCH --export ALL
#SBATCH -o worker${worker_num}.log

#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1

mkdir ${MYTMPDIR} -p
for i in \$(seq  \$SLURM_NTASKS ); do
        echo $i
        scripts/wq_worker_local.sh --cores 1 -s ${MYTMPDIR} --disk-threshold=0.002 --disk=3000 --memory-threshold=1000 -t 3600  -b 20 --memory=1000 $host:$port &
done
wait
EOF

sbatch $@ $cmd
rm $cmd
