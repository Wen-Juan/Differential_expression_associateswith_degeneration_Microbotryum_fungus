###This script is to build an index file for transcriptome, which is indeeded before running Kallisto mapping step.

module add UHTS/Analysis/kallisto/0.44.0

kallisto index -i Genome/transcriptome.idx Genome/transcriptome.fasta
