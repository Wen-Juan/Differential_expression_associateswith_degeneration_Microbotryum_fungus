cat a1_0k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a1_0k_TEnr.txt
cat a1_up5k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a1_up5k_TEnr.txt
cat a1_up10k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a1_up10k_TEnr.txt
cat a1_up15k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a1_up15k_TEnr.txt
cat a1_up20k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a1_up20k_TEnr.txt

wait

paste a1_0k_TEnr.txt a1_up5k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a1_up0-5k_TE_raw.txt
wait
cut -f1,5 a1_up0-5k_TE_raw.txt > a1_up0-5k_TE.txt

paste a1_up5k_TEnr.txt a1_up10k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a1_up5-10k_TE_raw.txt
wait
cut -f1,5 a1_up5-10k_TE_raw.txt > a1_up5-10k_TE.txt

paste a1_up10k_TEnr.txt a1_up15k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a1_up10-15k_TE_raw.txt
wait
cut -f1,5 a1_up10-15k_TE_raw.txt > a1_up10-15k_TE.txt

paste a1_up15k_TEnr.txt a1_up20k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a1_up15-20k_TE_raw.txt
wait
cut -f1,5 a1_up15-20k_TE_raw.txt > a1_up15-20k_TE.txt

wait

paste a1_0k_TEnr.txt a1_up0-5k_TE.txt a1_up5-10k_TE.txt a1_up10-15k_TE.txt a1_up15-20k_TE.txt > a1_upstream_allintervals.txt
