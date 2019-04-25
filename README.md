## Differential gene expression is associated with degeneration of mating-type chromosomes in the absence of sexual antagonism

In this project, we aim to investigate the association between differential gene expression and sequence degeneration on mating-type chromosome in an anther smut fungus species Microbotryum lychnidis-dioicae. The various forms of degenerative traits we are analyzing are non-synonymous substitution rate (dN), transposable element insertions (TEs), premature stop codons, intron content and GC content. 

To better navigate for the datasets and scripts, here are a brief README information. 

### Input folder including 8 sub-folders: 
1.1  /input/70percent_homologousgenes/
Within this folder, there are dataset files which were generated using various protein similarity thresholds: 30%, 70% and 85% using reciprocol best BLASTp Hit (RBBH). 

1.2  /input/Genomic_assignment/
Within this folder, there are genome annotation files in .gff3 formats, as well as gene locations and gene list of various evolutionary strata (young evolutionary strata include red and green strata; old evolutionary strata include purple, blue, orange and black strata).

1.3  /input/Kallisto_quantify_count/
Within this folder, there are count files (output file from Kallisto mapping), design and matrix files which were used for running downstream differential gene expression analysis using EdgeR. 

1.4  /input/dNdS/
Within this folder, there are files of dN, dS, and dN/dS data file within species M.lychnidis-dioicae.

1.5  /input/TE_degeneration/
Within this folder, there are files of TE insertions at various distance intervals both up and downstream 20kb, with and without oriented differences between given allele pairs based on the gene expression level (e.g. with higher expression or lower expression). There are also TE annotation files and TE consenseus sequence file (in Fasta format).

1.6  /input/Stopcodon/
Within this folder, there are files that show ratio of protein/coding sequence length between a1 and a2 alleles. There are also files for genes with unequal protein length between given alleles, and we then quantified the indel number, premature stop codon numbers and did comparisons across genomic compartments.

1.7  /input/Intron_content/
Within this folder, there are files that show both oriented and nonoriented intron length differences between DE and Non-DE genes. In addition, there are also coding sequence and protein fasta files for both a1 and a2 genomes, so we could use these files to generate the comparisons of gene prediction model bias against non-recombining regions.

1.8  /input/GCcontent/
Within this folder,there are files that show both oriended and nonoriended overall GC (GC0) and 3rd position (GC3) differences between given allele pairs.


### Script folder
Within this folder, there are scrips used for this projects. Including from RBBH.py to generate single copy homolog list between a1 and a2 haploid genomes, scripts to identify differential gene expression from Kallisto mapping to EdgeR analysis, prank.sh and PAML were used to calculate dN, dS etc data. Finally there are also R scripts for each degenerative traits for plotting some of the graphs. For each script, it always starts with a short explanation for what this script is for. Enjoy!
