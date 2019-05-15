# Within script folder, there are scripts used for this project. 


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

various R scripts: for each degenerative traits for plotting some of the graphs;scripts to identify differential gene expression from Kallisto pseudo-alignment mapping to EdgeR analysis; 

/script/GC_calculator_script/: scripts to calculate overall GC0 and GC3 percentage;

     bash GCcalculation.sh -a fasta_file 
     
....etc.



### For each script, it always starts with a short explanation for what this script is for.
