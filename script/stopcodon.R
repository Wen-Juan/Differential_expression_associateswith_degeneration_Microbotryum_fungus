#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load propotion data, modified codes on Jan.19.2019.
stopcodon_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/19jan2019/Mvsl_a1a2_exp_cds_protein_compart.txt', header = T)
str(stopcodon_ratio)

stopcodon_ratio_rmcentro <- subset(stopcodon_ratio,stopcodon_ratio$youngold != "Centro")
str(stopcodon_ratio_rmcentro)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protenlength_ratio_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=protratioa1a2, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,1.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of protein length (A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencelength_ratio_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=cdsratioa1a2, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,1.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length (A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencea1_devidprotein3times_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=cdsa1expest, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,2.2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencea2_devidprotein3times_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=cdsa2expest, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,2.2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load propotion data, modified codes on Jan.10.2019.
##protein length in 8 genomic compartments
stopcodon_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/stopcodon_de_ratio.txt', header = T)
str(stopcodon_ratio)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protenlength_ratio_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio, aes(factor(compartment), ratio, fill = bias)) + 
  scale_fill_manual(values = c("dodgerblue3","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

##protein length in 3 genomic compartments.
protlength_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/ratio_proteinlength.txt', header = T)
str(protlength_ratio)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protenlength_ratio_3PARTS.pdf", width=8, height=8)
ggplot(protlength_ratio, aes(factor(compartment), ratio, fill = bias)) + 
  scale_fill_manual(values = c("dodgerblue3","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

#cds length in 8 and 3 genomic compartments.
cds_eightparts <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/cds_8parts.txt', header = T)
str(cds_eightparts)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cdslength_ratio_8PARTS.pdf", width=8, height=8)
ggplot(cds_eightparts, aes(factor(compartment), ratio, fill = bias)) + 
  scale_fill_manual(values = c("dodgerblue3","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

cds_threeparts <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/cds_threeparts.txt', header = T)
str(cds_threeparts)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cdslength_ratio_3PARTS.pdf", width=8, height=8)
ggplot(cds_threeparts, aes(factor(compartment), ratio, fill = bias)) + 
  scale_fill_manual(values = c("dodgerblue3","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  ylim(0,1) +
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

#load the corresponding data files.
stopcodon <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/homolog_a1a2_cds_prot_length_removeTE_DEnonDEgeneexpree_intron_fi.txt', header = T)
str(stopcodon)

#########violin plot does not display message well####################
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cds_A1A2ratio_3PARTS.pdf", width=8, height=8)
ggplot(stopcodon, aes(chr1,cdsratio)) + 
  #geom_boxplot() +
  geom_violin(fill = "blue") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green","Red", "Orange","Black","Blue","Purple")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  labs(y='Coding sequence length ration (A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
#######################################################################

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_cdsa1a2ratio_3compartments.pdf", width=8, height=8)
ggplot(stopcodon, aes(x=chr1, y=cdsratio, fill=DE_status)) + 
  geom_boxplot(alpha=0.85) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("A2 biased","Not-biased","A1 biased"), name="Expression") + 
  ylim(0.5,1.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length (A1/A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_protiena1a2ratio_3compartments.pdf", width=8, height=8)
ggplot(stopcodon, aes(x=chr1, y=protratio, fill=DE_status)) + 
  geom_boxplot(alpha=0.85,position = position_dodge2(preserve = "single") ) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("A2 biased","Not-biased","A1 biased"), name="Expression") + 
  ylim(0.5,1.5) +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Ratio of protein length (A1/A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#scatter correlation
stopcodon_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/homolog_a1a2sep_cds_prot_length_removeTE_DEnonDEgeneexpree_tpm_fi.txt', header = T)
str(stopcodon_sep)

########not useful###############
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_8compartments_tpm.pdf", width=8, height=8)
ggplot(stopcodon_sep, aes(x=chr, y=log(meantpm), fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-3,7) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Log(meantpm)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
########not useful###############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_DEnonDE_protlength_3compartments_cor.pdf", width=15, height=8)
ggplot(stopcodon_sep, aes(x=chr1, y=protleng,fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"),labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1", "Not biased at A2","A1 biased at A1", "A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1600) +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Protein length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

y <- lm(log(meantpm+0.0001) ~ protleng*chr1-1, data = stopcodon_sep)
summary(y)

#############
lm(formula = log(meantpm + 1e-04) ~ protleng * chr1 - 1, data = stopcodon_sep)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.0232  -0.7605   0.0811   0.8646   7.5246 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
protleng          -9.713e-04  2.661e-05 -36.498   <2e-16 ***
  chr1aAutosome      2.812e+00  1.771e-02 158.802   <2e-16 ***
  chr1bPAR           2.644e+00  1.565e-01  16.900   <2e-16 ***
  chr1NRR            2.955e+00  1.183e-01  24.972   <2e-16 ***
  protleng:chr1bPAR  1.552e-04  2.666e-04   0.582    0.561    
protleng:chr1NRR   2.444e-04  1.670e-04   1.463    0.143    
#############