#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o samtools_bedfile_output.txt
#BSUB -e samtools_bedfile_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 20
#BSUB -R "span[ptile=20]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J samtools_bedfile.sh

module add UHTS/Analysis/samtools/1.8
 
samtools view -b -L M_inter_1389_onlygenes.bed ./Mint_output/accepted_hits.bam > accepted_hits_genes.bam
