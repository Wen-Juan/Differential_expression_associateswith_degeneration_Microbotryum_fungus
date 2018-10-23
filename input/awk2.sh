awk 'FNR==NR{seen[$4=1;next} ($1) in seen' A2_NRR_genes.txt A2_DEFC2_water_list.txt
