awk '$11 == 1 && $14 == "Auto" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonauto_cdsequal.txt
awk '$11 !=1 && $14 == "Auto" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonauto_cdsnotequal.txt
awk '$12 == 1 && $14 == "Auto" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonauto_proteinequal.txt
awk '$12 != 1 && $14 == "Auto" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonauto_proteinnotequal.txt

awk '$11 == 1 && $14 == "bPAR" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonbPAR_cdsequal.txt
awk '$11 !=1 && $14 == "bPAR" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonbPAR_cdsnotequal.txt
awk '$12 == 1 && $14 == "bPAR" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonbPAR_proteinequal.txt
awk '$12 != 1 && $14 == "bPAR" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonbPAR_proteinnotequal.txt

awk '$11 == 1 && $14 == "ColorStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonColorStrata_cdsequal.txt
awk '$11 !=1 && $14 == "ColorStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonColorStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "ColorStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonColorStrata_proteinequal.txt
awk '$12 != 1 && $14 == "ColorStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonColorStrata_proteinnotequal.txt

awk '$11 == 1 && $14 == "OldStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonOldStrata_cdsequal.txt
awk '$11 !=1 && $14 == "OldStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonOldStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "OldStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonOldStrata_proteinequal.txt
awk '$12 != 1 && $14 == "OldStrata" && $25 == "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > nonOldStrata_proteinnotequal.txt



awk '$11 == 1 && $14 == "Auto" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEauto_cdsequal.txt
awk '$11 !=1 && $14 == "Auto" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEauto_cdsnotequal.txt
awk '$12 == 1 && $14 == "Auto" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEauto_proteinequal.txt
awk '$12 != 1 && $14 == "Auto" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEauto_proteinnotequal.txt

awk '$11 == 1 && $14 == "bPAR" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEbPAR_cdsequal.txt
awk '$11 !=1 && $14 == "bPAR" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEbPAR_cdsnotequal.txt
awk '$12 == 1 && $14 == "bPAR" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEbPAR_proteinequal.txt
awk '$12 != 1 && $14 == "bPAR" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEbPAR_proteinnotequal.txt

awk '$11 == 1 && $14 == "ColorStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEColorStrata_cdsequal.txt
awk '$11 !=1 && $14 == "ColorStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEColorStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "ColorStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEColorStrata_proteinequal.txt
awk '$12 != 1 && $14 == "ColorStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEColorStrata_proteinnotequal.txt

awk '$11 == 1 && $14 == "OldStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEOldStrata_cdsequal.txt
awk '$11 !=1 && $14 == "OldStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEOldStrata_cdsnotequal.txt
awk '$12 == 1 && $14 == "OldStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEOldStrata_proteinequal.txt
awk '$12 != 1 && $14 == "OldStrata" && $25 != "NON"' Mvsl_a1a2_exp_cds_protein_compart.txt | wc -l > DEOldStrata_proteinnotequal.txt
