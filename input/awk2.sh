awk 'FNR==NR{seen[$4]=1;next} ($1) in seen' A1_NRR_genes.txt UP_A1_FDR005_genes.txt
