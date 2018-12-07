#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
stopcodon <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/homolog_a1a2_cds_prot_length_removeTE_DEnonDEgeneexpree_intron_fi.txt', header = T)
str(stopcodon)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1A2_cdsratio_a1a2_3PARTS.pdf", width=8, height=8)
ggplot(stopcodon, aes(chr1,cdsratio)) + 
  #geom_boxplot() +
  geom_violin(fill = "blue") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green","Red", "Orange","Black","Blue","Purple")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  labs(y='Coding sequence length ration (A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(stopcodon, aes(x=chr1, y=cdsratio, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("Down","Up", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  #scale_x_discrete(labels=c("up:10-22k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_cdsratio_3compartments.pdf", width=8, height=8)
ggplot(stopcodon, aes(x=logFC.A1.A2, y=cdsratio, fill=DE_status)) + 
  geom_boxplot(alpha=0.85,position = position_dodge2(preserve = "single") ) +
  facet_grid(~chr1) +
  #ylim(0.9,1.1) +
  labs(y='Ration of protein length (A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

ggplot(stopcodon, aes(x=logFC.A1.A2, y=protratio, fill=DE_status)) + 
 # scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("Down","Up", "not bias"), name="DE expression") + 
  geom_boxplot(alpha=0.85,position = position_dodge2(preserve = "single") ) +
  facet_grid(~chr1) +
  ylim(0,0.5)
labs(y='Protein length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

#scatter correlation
stopcodon_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/homolog_a1a2sep_cds_prot_length_removeTE_DEnonDEgeneexpree_tpm_fi.txt', header = T)
str(stopcodon_sep)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_8compartments_tpm.pdf", width=8, height=8)
ggplot(stopcodon_sep, aes(x=chr, y=log(meantpm), fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-3,7) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Log(meantpm)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_DEnonDE_protlength_3compartments_cor.pdf", width=15, height=8)
ggplot(stopcodon_sep, aes(x=log(meantpm+0.00001), y=protleng,fill=interaction(haploid,DE_status))) + 
  geom_boxplot(alpha=0.85, position = position_dodge2(),outlier.shape=NA,width=0.2) +
  facet_grid(~chr1) +
  ylim(0,1600) +
 # xlim(0,100) +
  labs(y='Protein length') +
  labs(x='Log(meanTPM)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2sep_DEnonDE_protlength_3compartment_scatter1.pdf", width=8, height=8)
ggplot(stopcodon_sep, aes(x=log(meantpm), y=protleng,color=interaction(haploid,DE_status))) + 
  geom_point (alpha=0.5)+
  facet_grid(~chr1) +
  #geom_smooth(method=lm,se=FALSE) +
  ylim(0,4000) +
  labs(y='Protein length') +
  labs(x='Log(meanTPM)') +
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