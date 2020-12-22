#!/bin/bash
date
module list

REDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit
EDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit_edlib
REF=Chromosome_2890043_3890042_0/Chromosome_2890043_3890042_0.fasta

cat Chromosome_2890043_3890042_0/Chromosome_2890043_3890042_0.fasta Chromosome_2890043_3890042_0/mutated_97_perc.fasta Chromosome_2890043_3890042_0/mutated_94_perc.fasta Chromosome_2890043_3890042_0/mutated_90_perc.fasta Chromosome_2890043_3890042_0/mutated_80_perc.fasta Chromosome_2890043_3890042_0/mutated_70_perc.fasta Chromosome_2890043_3890042_0/mutated_60_perc.fasta > query.fa 

QUERY=query.fa

/usr/bin/time $EDIT -q $QUERY -t $REF -l -1 -m g -a MEM &> log_edit &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time $REDIT -q $QUERY -t $REF -l 20 -m g -a MEM &> log_redit20 &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time $REDIT -q $QUERY -t $REF -l 10 -m g -a MEM &> log_redit10 &

wait
date
echo "DONE!"
