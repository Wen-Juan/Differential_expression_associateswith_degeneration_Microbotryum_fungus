 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_2kupstream_TE_number.txt &
 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_10kupstream__TE_number.txt
 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_20kupstream__TE_number.txt

 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_2downpstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_2kupstream_TE_number.txt &
 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_10kupstream__TE_number.txt
 grep -f Mvsl_DEnonDE_A1_genelist.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A1_20kupstream__TE_number.txt

 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_2kupstream_TE_number.txt &
 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_10kupstream__TE_number.txt
 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_20kupstream__TE_number.txt

 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_2downpstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_2kupstream_TE_number.txt &
 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_10kupstream_TE_number.txt
 grep -f Mvsl_DEnonDE_A2_genelist.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 |less | sort | uniq -c | sort -n | tr ':' '\t' | tr ';' '\t' | sed -e $'s/ /\t/g' | cut -f4,6 > Mvsl_DEnonDE_A2_20kupstream__TE_number.txt
