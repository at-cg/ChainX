#!/bin/bash
date
module list

REDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit
EDIT=/global/homes/c/cjain7/work/projects/redit/code/rEDIT/redit_edlib

REF=E_coli_DH1/e_coli_DH1.fasta
QUERY=E_coli_DH1/mason_illumina_reads/500bp/mutated_97_perc.fasta

/usr/bin/time $REDIT -q $QUERY -t $REF -l 10 -m g -a MEM --naive &> log_dp-g10 & 
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time $REDIT -q $QUERY -t $REF -l 10 -m g -a MEM &> log_g10 &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time $REDIT -q $QUERY -t $REF -l 10 -m sg -a MEM --naive &> log_dp-sg10 & 
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time $REDIT -q $QUERY -t $REF -l 10 -m sg -a MEM &> log_sg10 & 

wait
date
echo "DONE!"
