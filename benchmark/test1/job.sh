#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --constraint=haswell
#SBATCH --output=BATCH_OUTPUT
#SBATCH --error=BATCH_OUTPUT

date
module list

REDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit
EDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit_edlib

REF=Chromosome_2890043_3890042_0/Chromosome_2890043_3890042_0.fasta

cat Chromosome_2890043_3890042_0/Chromosome_2890043_3890042_0.fasta Chromosome_2890043_3890042_0/mutated_97_perc.fasta Chromosome_2890043_3890042_0/mutated_94_perc.fasta Chromosome_2890043_3890042_0/mutated_90_perc.fasta Chromosome_2890043_3890042_0/mutated_80_perc.fasta Chromosome_2890043_3890042_0/mutated_70_perc.fasta Chromosome_2890043_3890042_0/mutated_60_perc.fasta > query.fa 

QUERY=query.fa

/usr/bin/time srun -N 1 $EDIT -q $QUERY -t $REF -l -1 -m g &> log_edit &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time srun -N 1 $REDIT -q $QUERY -t $REF -l 20 -m g &> log_redit20 &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time srun -N 1 $REDIT -q $QUERY -t $REF -l 10 -m g &> log_redit10 &

wait
date
echo "DONE!"
