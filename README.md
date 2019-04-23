## Differential gene expression is associated with degeneration of mating-type chromosomes in the absence of sexual antagonism

In this project, we aim to investigate the association between differential gene expression and sequence degeneration on mating-type chromosome in an anther smut fungus species Microbotryum lychnidis-dioicae. The various forms of degenerative traits we are analyzing are synonymous substitution rate (dN), transposable element insertions (TEs), premature stop codons, intron content and GC content. 

To better navigate for the datasets and scripts, here are a brief README information. 

There are in total 8 different sub-folders for this project input folder: 
1.1/input/70percent_homologousgenes/
Within this folder, there are dataset files which were generated using various protein similarity thresholds: 30%, 70% and 85% using reciprocol best BLASTp Hit (RBBH). 

1.2/input/Genomic_assignment/
Within this folder, there are genome annotation files in .gff3 formats, as well as genome locations and gene list of various evolutionary strata (young evolutionary strata include red and green strata; old evolutionary strata include purple, blue, orange and black strata).

1.3/input/Kallisto_quantify_count/
Within this folder there are count file (output file frim Kallisto mapping), design and matrix files which were used for running downstream differential gene expression analysis using EdgeR. 

1.4/input/dNdS/
Within this folder there are files of dN, dS, and dN/dS data file.

1.5./input/TE_degeneration/
Within this folder there are files of TE insertions at various distance intervals both up and downstream 20kb, with and without oriented differences between given allele pairs. There are also TE annotation files and TE consenseus sequence file (in Fasta format).

1.6/input/Genomic_assignment/

1.7/input/Genomic_assignment/

1.8/input/Genomic_assignment/
