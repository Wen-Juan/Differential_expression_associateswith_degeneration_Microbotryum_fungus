# Within script folder, there are scripts used for this project. 

### For each script, it always starts with a short explanation for what this script is for.

## Scripts include:

trimmomatic.sh: trim raw reads to remove low quality reads, adapters etc. using default parameters.

    bash trimmomatic.sh
fastqc.sh: assess read quality before and after trimming.

    bash fastqc.sh
RBBH.py: generate single-copy homolog list between a1 and a2 haploid genomes. In the instruction of this python code, you can define the threshold of similarity between homologous protein pairs, coverage etc.

    python RBBH.py

Extract_paired_sequence.py, or extract_orthologs_seqs.py: scripts to extract paired sequences from a text file of two lists and two fasta files.
    
    python Extract_paired_sequence.py
    python extract_orthologs_seqs.py
extract_seq_from_fasta.py: scripts to extract certain sequences with indicated headers.
    
    python extract_seq_from_fasta.py
prank.sh: to align sequence in a pair-wise fashion for input of PAML.

    bash prank.sh

loop_codeml_se.sh: to calculate dN, dS, dN/dS gene evolutionary rate;

    bash loop_codeml_se.sh
kallisto_index.sh: build index of transcriptome for kallisto mapping.

    bash kallisto_count.sh
kallisto_count.sh: kallisto to map RNAseq reads against transcriptome.

    bash kallisto_count.sh
kallisto_matrix.sh: kallisto matrix script.

    bash kallisto_matrix.sh
extract_intron_size_from_gff_files.pl: script to calculate intron number and length from gff file.

    perl extract_intron_size_from_gff_files.pl
/script/GC_calculator_script/: scripts to calculate overall GC0 and GC3 percentage;

     bash GCcalculation.sh -a fasta_file 
various R scripts, edgeR.r: scripts to identify differential gene expression from Kallisto pseudo-alignment mapping to EdgeR analysis. 
     
     edgeR.r

below R scripts are for each degenerative traits, e.g. dN, dS, TEs, premature stop codons, GC and introns for plotting the graphs and analysis.
     
     dNdS.R
     TE_degeneration.R
     Premature_stopcodon.R
     GC_content.R
     intron_content.R
     
