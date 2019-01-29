awk '$8 == "Auto" && $20 == "DE"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_auto_DE.txt
awk '$8 == "Auto" && $20 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_auto_NON.txt

awk '$8 == "bPAR" && $20 == "DE"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_par_DE.txt
awk '$8 == "bPAR" && $20 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_par_NON.txt

awk '$8 == "ColorStrata" && $20 == "DE"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_colorstrata_DE.txt
awk '$8 == "ColorStrata" && $20 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_colorstrata_NON.txt

awk '$8 == "OldStrata" && $20 == "DE"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_oldstrata_DE.txt
awk '$8 == "OldStrata" && $20 == "NON"' Mvsl_dnds_exp_gencomp.txt | wc -l > samplesize_oldstrata_NON.txt
