awk '$11 == 1 && $14 == "Auto"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > auto_cdsequal.txt
awk '$11 !=1 && $14 == "Auto"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > auto_cdsnotequal.txt
awk '$12 == 1 && $14 == "Auto"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > auto_proteinequal.txt
awk '$12 != 1 && $14 == "Auto"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > auto_proteinnotequal.txt

awk '$11 == 1 && $14 == "bPAR"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > bPAR_cdsequal.txt
awk '$11 !=1 && $14 == "bPAR"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > bPAR_cdsnotequal.txt
awk '$12 == 1 && $14 == "bPAR"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > bPAR_proteinequal.txt
awk '$12 != 1 && $14 == "bPAR"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > bPAR_proteinnotequal.txt

awk '$11 == 1 && $14 == "ColorStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > ColorStrata_cdsequal.txt
awk '$11 !=1 && $14 == "ColorStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > ColorStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "ColorStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > ColorStrata_proteinequal.txt
awk '$12 != 1 && $14 == "ColorStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > ColorStrata_proteinnotequal.txt

awk '$11 == 1 && $14 == "OldStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > OldStrata_cdsequal.txt
awk '$11 !=1 && $14 == "OldStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > OldStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "OldStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > OldStrata_proteinequal.txt
awk '$12 != 1 && $14 == "OldStrata"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > OldStrata_proteinnotequal.txt
