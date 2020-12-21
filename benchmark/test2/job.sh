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

REF=E_coli_DH1/e_coli_DH1.fasta

cat E_coli_DH1/mason_illumina_reads/10kbp/e_coli_DH1_illumina_1x10000.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_97_perc.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_94_perc.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_90_perc.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_80_perc.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_70_perc.fasta E_coli_DH1/mason_illumina_reads/10kbp/mutated_60_perc.fasta > query.fa

QUERY=query.fa

/usr/bin/time srun -N 1 $EDIT -q $QUERY -t $REF -l -1 -m sg &> log_edit &
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time srun -N 1 $REDIT -q $QUERY -t $REF -l 20 -m sg &> log_redit20 & 
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
/usr/bin/time srun -N 1 $REDIT -q $QUERY -t $REF -l 10 -m sg &> log_redit10 &

wait
date
echo "DONE!"
