# Within script folder, there are scripts used for this project. 



## Scripts include:

trimmomatic.sh: trim raw reads to remove low quality reads, adapters etc. using default parameters.

    bash trimmomatic.sh

RBBH.py: generate single-copy homolog list between a1 and a2 haploid genomes, scripts to identify differential gene expression from Kallisto pseudo-alignment mapping to EdgeR analysis; In the instruction of this python code, you can define the threshold of similarity between homologous protein pairs, coverage etc.

    python RBBH.py

prank.sh: to align sequence in a pair-wise fashion for input of PAML.

    bash prank.sh

loop_codeml_se.sh: to calculate dN, dS, dN/dS gene evolutionary rate;

    bash loop_codeml_se.sh

various R scripts: for each degenerative traits for plotting some of the graphs;

/script/GC_calculator_script/: scripts to calculate overall GC0 and GC3 percentage;

     python 
     
....etc.



### For each script, it always starts with a short explanation for what this script is for.
