cat a2_0k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a2_0k_TEnr.txt
cat a2_down5k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a2_down5k_TEnr.txt
cat a2_down10k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a2_down10k_TEnr.txt
cat a2_down15k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a2_down15k_TEnr.txt
cat a2_down20k_output.txt | sed 's/Name=//g' | cut -f9,10 | tr ';' \\t | cut -f2,4 | sort -k1,1 > a2_down20k_TEnr.txt

wait

paste a2_0k_TEnr.txt a2_down5k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a2_down0-5k_TE_raw.txt
wait
cut -f1,5 a2_down0-5k_TE_raw.txt > a2_down0-5k_TE.txt

paste a2_down5k_TEnr.txt a2_down10k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a2_down5-10k_TE_raw.txt
wait
cut -f1,5 a2_down5-10k_TE_raw.txt > a2_down5-10k_TE.txt

paste a2_down10k_TEnr.txt a2_down15k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a2_down10-15k_TE_raw.txt
wait
cut -f1,5 a2_down10-15k_TE_raw.txt > a2_down10-15k_TE.txt

paste a2_down15k_TEnr.txt a2_down20k_TEnr.txt | awk 'NR == 1 { $5 = "diff" } NR >= 1 { $5 = $4 - $2 } 1' | tr ' ' \\t > a2_down15-20k_TE_raw.txt
wait
cut -f1,5 a2_down15-20k_TE_raw.txt > a2_down15-20k_TE.txt

wait

paste a2_down0-5k_TE.txt a2_down5-10k_TE.txt a2_down10-15k_TE.txt a2_down15-20k_TE.txt > a2_downstream_allintervals.txt
