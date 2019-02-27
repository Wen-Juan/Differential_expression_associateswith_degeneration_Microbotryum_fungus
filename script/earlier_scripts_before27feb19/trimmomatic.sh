#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o trim_output.txt
#BSUB -e trim_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 15
#BSUB -R "span[ptile=15]"
#BSUB -R "rusage[mem=40000]"
#BSUB -M 40000000
#BSUB -J trimmomatic.sh

module add UHTS/Analysis/trimmomatic/0.36

ADAPTERS="/scratch/beegfs/weekly/wjma/RNAseq_Microbotryum/TruSeq2-PE.fa"

for f in *_1.fq.gz

do
     trimmomatic PE -phred33 \
     ${f%%_1.fq.gz}_1.fq.gz \
     ${f%%_1.fq.gz}_2.fq.gz \
     ${f%%_1.fq.gz}_pairedR1.fq.gz \
     ${f%%_1.fq.gz}_unpairedR1.fq.gz \
     ${f%%_1.fq.gz}_pairedR2.fq.gz \
     ${f%%_1.fq.gz}_unpairedR2.fq.gz \
     ILLUMINACLIP:$ADAPTERS:2:30:10 \
     LEADING:3 \
     TRAILING:3 \
     SLIDINGWINDOW:4:15 \
     MINLEN:36 &> ${f%%1combined.fq.gz}.trim.log

done
