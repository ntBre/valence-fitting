#!/bin/bash -x

USERNAME=$(whoami)
MYTMPDIR="/tmp/${USERNAME}"

mkdir -p $MYTMPDIR
cd "/tmp/${USERNAME}"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

work_queue_worker $@
