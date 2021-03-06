#This script is to build matrix of multiple RNAseq samples with read counts, using output from Kallisto for instance.

##before launching the matrix perl script, you need to load edgeR and qvalue packages from R.
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")

module add UHTS/Assembler/trinityrnaseq/2.5.1
module add UHTS/Analysis/kallisto/0.44.0
module add R/3.5.1

/software/UHTS/Assembler/trinityrnaseq/2.5.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  \
--gene_trans_map none \
--out_prefix Haploid_rich_ \
--name_sample_by_basedir \
/A1/SRR3624826/abundance.tsv \
/A1/SRR3710483/abundance.tsv \
/A2/SRR3710080/abundance.tsv \
/A2/SRR3710485/abundance.tsv
