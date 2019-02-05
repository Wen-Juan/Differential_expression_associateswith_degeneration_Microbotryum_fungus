awk '$3 == "Auto" && $20 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > auto_cdslength_DE_not1.txt
awk '$3 == "Auto" && $20 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > auto_cdslength_NON_not1.txt

awk '$3 == "bPAR" && $20 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > bPAR_cdslength_DE_not1.txt
awk '$3 == "bPAR" && $20 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > bPAR_cdslength_NON_not1.txt

awk '$3 == "ColorStrata" && $20 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > ColorStrata_cdslength_DE_not1.txt
awk '$3 == "ColorStrata" && $20 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > ColorStrata_cdslength_NON_not1.txt

awk '$3 == "OldStrata" && $20 != 1 && $17 !="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > OldStrata_cdslength_DE_not1.txt
awk '$3 == "OldStrata" && $20 != 1 && $17 =="NON"' Mvsl_a1a2_exp_gencompt_protlength_fi.txt > OldStrata_cdslength_NON_not1.txt
