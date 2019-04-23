# LSBATCH: User input
#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o fastqc2_output.txt
#BSUB -e fastqc2_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 20
#BSUB -R "span[ptile=20]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J fastqc2.sh

module add UHTS/Quality_control/fastqc/0.11.7

fastqc -o ./raw_qc/ M1_M1_1.fq.gz M1_M1_2.fq.gz &
fastqc -o ./raw_qc/ M1_M2_1.fq.gz M1_M2_2.fq.gz &
fastqc -o ./raw_qc/ M1_M3_1.fq.gz M1_M3_2.fq.gz
