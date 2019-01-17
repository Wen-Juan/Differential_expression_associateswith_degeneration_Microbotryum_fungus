###check total DE and non-DE gene number with TE insertion sites.
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_auto_a1nr.txt &
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_auto_a1nr.txt
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_auto_a1nr.txt
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_auto_a1nr.txt &
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_auto_a1nr.txt
grep -f Mvsl_a1bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_auto_a1nr.txt

grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_auto_a1nr.txt &
grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_auto_a1nr.txt
grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_auto_a1nr.txt
grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_auto_a1nr.txt &
grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_auto_a1nr.txt
grep -f Mvsl_a2bias_autosomegenes.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_auto_a1nr.txt

grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kupstream_auto_2nr.txt &
grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kupstream_auto_a2nr.txt
grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kupstream_auto_a2nr.txt
grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_2kdownstream_auto_a2nr.txt &
grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_10kdownstream_auto_a2nr.txt
grep -f Mvsl_a1bias_autosomeA2.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a1bias_20kdownstream_auto_a2nr.txt

grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kupstream_auto_a2nr.txt &
grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kupstream_auto_a2nr.txt
grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kupstream_auto_a2nr.txt
grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_2kdownstream_auto_a2nr.txt &
grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_10kdownstream_auto_a2nr.txt
grep -f Mvsl_a2bias_autosomesA2.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_a2bias_20kdownstream_auto_a2nr.txt

grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_auto_a1nr.txt &
grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_auto_a1nr.txt
grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_auto_a1nr.txt
grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_auto_a1nr.txt &
grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_auto_a1nr.txt
grep -f Mvsl_nonDE_auto_A1genes_list.txt MvSl_A1_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_auto_a1nr.txt

grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_2kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kupstream_auto_a2nr.txt &
grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_10kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kupstream_auto_a2nr.txt
grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_20kupstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kupstream_auto_a2nr.txt
grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_2kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_2kdownstream_auto_a2nr.txt &
grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_10kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_10kdownstream_auto_a2nr.txt
grep -f Mvsl_nonDE_auto_A2genes_list.txt MvSl_A2_20perc_overlap_TE_20kdownstream_TEs.txt | cut -f9 | uniq | wc -l > Mvsl_nonDE_20kdownstream_auto_a2nr.txt
