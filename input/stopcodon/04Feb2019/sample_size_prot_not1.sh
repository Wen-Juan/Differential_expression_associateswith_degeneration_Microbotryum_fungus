awk '$3 == "Auto" && $16 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > auto_protlength_DE_not1.txt
awk '$3 == "Auto" && $16 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > auto_protlength_NON_not1.txt

awk '$3 == "bPAR" && $16 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > bPAR_protlength_DE_not1.txt
awk '$3 == "bPAR" && $16 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > bPAR_protlength_NON_not1.txt

awk '$3 == "ColorStrata" && $16 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > ColorStrata_protlength_DE_not1.txt
awk '$3 == "ColorStrata" && $16 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > ColorStrata_protlength_NON_not1.txt

awk '$3 == "OldStrata" && $16 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > OldStrata_protlength_DE_not1.txt
awk '$3 == "OldStrata" && $16 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > OldStrata_protlength_NON_not1.txt
