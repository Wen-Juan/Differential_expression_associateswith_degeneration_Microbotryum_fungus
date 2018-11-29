###check total DE and non-DE gene number with TE insertion sites.
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_a1nr.txt &
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_a1nr.txt
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_a1nr.txt
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_a1nr.txt &
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_a1nr.txt
grep -f Mvsl_a1bias_gene_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_a1nr.txt

grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_a1nr.txt &
grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_a1nr.txt
grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_a1nr.txt
grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_a1nr.txt &
grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_a1nr.txt
grep -f Mvsl_a2bias_gene_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_a1nr.txt

grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_a2nr.txt &
grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_a2nr.txt
grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_a2nr.txt
grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_a2nr.txt &
grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_a2nr.txt
grep -f MvSl_DE_a1bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_a2nr.txt

grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_a2nr.txt &
grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_a2nr.txt
grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_a2nr.txt
grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_a2nr.txt &
grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_a2nr.txt
grep -f MvSl_DE_a2bias_a2genelist.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_a2nr.txt

grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_a1nr.txt &
grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_a1nr.txt
grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_a1nr.txt
grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_a1nr.txt &
grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_a1nr.txt
grep -f Mvsl_nonDE_A1list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_a1nr.txt

grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_a2nr.txt &
grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_a2nr.txt
grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_a2nr.txt
grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_a2nr.txt &
grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_a2nr.txt
grep -f Mvsl_A2_nonDE_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_a2nr.txt
