#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o tophat_mint_output.txt
#BSUB -e tophat_mint_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 20
#BSUB -R "span[ptile=20]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J tophat_mint.sh

module add UHTS/Aligner/tophat/2.1.1

module add UHTS/Aligner/bowtie2/2.3.4.1

module add UHTS/Analysis/samtools/1.8

tophat -p 20 -G M_inter_1389.gff3 -o ./Mint_output M_inter_1389 M1_M1_R1.fq.gz,M1_M1_R2.fq.gz
