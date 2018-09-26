#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o Mint_idxstats_output.txt
#BSUB -e Mint_idxstats_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 20
#BSUB -R "span[ptile=20]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J Mint_idxstats.sh

module add UHTS/Analysis/samtools/1.8

samtools sort -o Mint_accepted_hits_genes_sorted.bam Mint_accepted_hits_genes.bam

wait

samtools index Mint_accepted_hits_genes_sorted.bam

wait

samtools idxstats Mint_accepted_hits_genes_sorted.bam
