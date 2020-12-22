#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --constraint=haswell
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

srun -N 1 ./job.sh
