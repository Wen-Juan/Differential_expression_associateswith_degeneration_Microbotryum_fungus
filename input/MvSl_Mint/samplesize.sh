awk '$3 == "Auto" && $8 =="Vhighmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_Auto_Vhighmutations.txt
awk '$3 == "Auto" && $8 =="Lowmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_Auto_Lowmutations.txt
awk '$3 == "Auto" && $8 =="Neutral"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_Auto_Neutral.txt

awk '$3 == "bPAR" && $8 =="Vhighmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_bPAR_Vhighmutations.txt
awk '$3 == "bPAR" && $8 =="Lowmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_bPAR_Lowmutations.txt
awk '$3 == "bPAR" && $8 =="Neutral"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_bPAR_Neutral.txt

awk '$3 == "ColorStrata" && $8 =="Vhighmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_ColorStrata_Vhighmutations.txt
awk '$3 == "ColorStrata" && $8 =="Lowmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_ColorStrata_Lowmutations.txt
awk '$3 == "ColorStrata" && $8 =="Neutral"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_ColorStrata_Neutral.txt

awk '$3 == "OldStrata" && $8 =="Vhighmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_OldStrata_Vhighmutations.txt
awk '$3 == "OldStrata" && $8 =="Lowmutations"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_OldStrata_Lowmutations.txt
awk '$3 == "OldStrata" && $8 =="Neutral"' Mvsl_Mvsv_dnds_exp_all_fi_sep.txt | wc -l > samplesize_OldStrata_Neutral.txt
