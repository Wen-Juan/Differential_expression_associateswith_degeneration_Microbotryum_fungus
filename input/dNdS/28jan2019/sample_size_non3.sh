awk '$8 == "Auto" && $19 == "Down"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_auto_Down.txt
awk '$8 == "Auto" && $19 == "Up"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_auto_Up.txt
awk '$8 == "Auto" && $19 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_auto_NON3.txt

awk '$8 == "bPAR" && $19 == "Down"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_bPAR_Down.txt
awk '$8 == "bPAR" && $19 == "Up"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_bPAR_Up.txt
awk '$8 == "bPAR" && $19 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_bPAR_NON3.txt

awk '$8 == "ColorStrata" && $19 == "Down"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_ColorStrata_Down.txt
awk '$8 == "ColorStrata" && $19 == "Up"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_ColorStrata_Up.txt
awk '$8 == "ColorStrata" && $19 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_ColorStrata_NON3.txt

awk '$8 == "OldStrata" && $19 == "Down"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_OldStrata_Down.txt
awk '$8 == "OldStrata" && $19 == "Up"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_OldStrata_Up.txt
awk '$8 == "OldStrata" && $19 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_OldStrata_NON3.txt
