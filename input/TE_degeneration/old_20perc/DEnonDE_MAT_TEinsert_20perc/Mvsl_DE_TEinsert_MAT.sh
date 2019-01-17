###check total DE and non-DE gene number with TE insertion sites.
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_MAT_a1nr.txt &
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_MAT_a1nr.txt
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_MAT_a1nr.txt
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_MAT_a1nr.txt &
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_MAT_a1nr.txt
grep -f Mvsl_a1bias_MATgenes.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_MAT_a1nr.txt

grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_MAT_a1nr.txt &
grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_MAT_a1nr.txt
grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_MAT_a1nr.txt
grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_MAT_a1nr.txt &
grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_MAT_a1nr.txt
grep -f Mvsl_a2bias_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_MAT_a1nr.txt

grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_MAT_2nr.txt &
grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_MAT_a2nr.txt
grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_MAT_a2nr.txt
grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_MAT_a2nr.txt &
grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_MAT_a2nr.txt
grep -f Mvsl_a1bias_MATgenesA2.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_MAT_a2nr.txt

grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_MAT_a2nr.txt &
grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_MAT_a2nr.txt
grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_MAT_a2nr.txt
grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_MAT_a2nr.txt &
grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_MAT_a2nr.txt
grep -f Mvsl_a2bias_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_MAT_a2nr.txt

grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_MAT_a1nr.txt &
grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_MAT_a1nr.txt
grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_MAT_a1nr.txt
grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_MAT_a1nr.txt &
grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_MAT_a1nr.txt
grep -f Mvsl_nonDE_MAT_a1genes_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_MAT_a1nr.txt

grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_MAT_a2nr.txt &
grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_MAT_a2nr.txt
grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_MAT_a2nr.txt
grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_MAT_a2nr.txt &
grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_MAT_a2nr.txt
grep -f Mvsl_nonDE_MAT_a2genes_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_MAT_a2nr.txt
