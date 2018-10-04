# LSBATCH: User input
#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o kallisto_mint_output.txt
#BSUB -e kallisto_mint_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 20
#BSUB -R "span[ptile=20]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J kallisto_mint.sh

module add UHTS/Analysis/kallisto/0.43.0

for f in *_pairedR1.fq.gz
do
kallisto quant -i /scratch/beegfs/weekly/wjma/RNAseq_Microbotryum/genome_gene/Mint/M_inter_1389.idx -o ./output/${f%%_pairedR1.fq.gz} -b 1000 ${f%%_pairedR1.fq.gz}_pairedR1.fq.gz ${f%%_pairedR1.fq.gz}_pairedR2.fq.gz
done
