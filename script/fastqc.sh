#This script is to assess quality control of raw or trimmed paired-end fastq read files.

module add UHTS/Quality_control/fastqc/0.11.7

fastqc -o ./raw_qc/ M1_M1_1.fq.gz M1_M1_2.fq.gz &
fastqc -o ./raw_qc/ M1_M2_1.fq.gz M1_M2_2.fq.gz &
fastqc -o ./raw_qc/ M1_M3_1.fq.gz M1_M3_2.fq.gz
