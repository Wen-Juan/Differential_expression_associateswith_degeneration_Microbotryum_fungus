###check total DE and non-DE gene number with TE insertion sites.
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_MAT_a1nr.txt &
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_richmedium_a1nr.txt &
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a1bias_genes_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_richmedium_a1nr.txt

grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_richmedium_a1nr.txt &
grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_richmedium_a1nr.txt &
grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_richmedium_a1nr.txt
grep -f Mvsl_richmedium_a2bias_genes_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_richmedium_a1nr.txt

grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_richmedium_2nr.txt &
grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_richmedium_a2nr.txt
grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_richmedium_a2nr.txt
grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_richmedium_a2nr.txt &
grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_richmedium_a2nr.txt
grep -f Mvsl_richmedium_a1bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20Mvsl_a2bias_20kdownstream_MAT_a2nr_richmedium_a2nr.txt

grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_richmedium_a2nr.txt &
grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_richmedium_a2nr.txt
grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_richmedium_a2nr.txt
grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_richmedium_a2nr.txt &
grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_richmedium_a2nr.txt
grep -f Mvsl_richemedium_a2bias_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_MAT_a2nr.txt

grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_MAT_a1nr.txt &
grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_MAT_a1nr.txt
grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_MAT_a1nr.txt
grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_MAT_a1nr.txt &
grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_MAT_a1nr.txt
grep -f Mvsl_richmedium_nonDE_genesa1_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_MAT_a1nr.txt

grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_MAT_a2nr.txt &
grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_MAT_a2nr.txt
grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_MAT_a2nr.txt
grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_MAT_a2nr.txt &
grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_MAT_a2nr.txt
grep -f Mvsl_richmedium_nonDE_genesa2_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_MAT_a2nr.txt
