#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load data on 13feb2019 for directional ratio, and randomdize non-DE genes directions for ratio calculation.
non_DE_protein <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/Mvsl_a1a2_proteinlength_nonDE.txt', header = T)
str(non_DE_protein)

prot_ratio <- data.frame(non_DE_protein$a1prot, non_DE_protein$a2prot)
cds_ratio <- data.frame(non_DE_protein$a1cds, non_DE_protein$a2cds)

prot_ratio_rand <- randomizeMatrix(prot_ratio,null.model = "richness",iterations = 1000)
cds_ratio_rand <- randomizeMatrix(cds_ratio,null.model = "richness",iterations = 1000)

total_data <- cbind(non_DE_protein, prot_ratio_rand,cds_ratio_rand)
write.table(total_data,file = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/Mvsl_a1a2_nonDE_prot_randomdized.txt",quote=F, row.names=T, sep='\t')

diff_prot_length_rand <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/Mvsl_a1a2_exp_gencompt_protlength_fi_randomdize.txt', header = T)
str(diff_prot_length_rand)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protenlength_ratio_rand_youngold.pdf", width=8, height=8)
ggplot(diff_prot_length_rand, aes(x=ratioprot, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(0,2) +
  labs(x='Ratio of protein length', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m1 <- lm(abs ~ DE2/ratioprot-1, data = diff_prot_length_rand)
summary(y_m1)
###################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            -1.13821    0.25928  -4.390 1.15e-05 ***
  DE2NON            0.26189    0.26277   0.997    0.319    
DE2DE:ratioprot   3.02000    0.25820  11.696  < 2e-16 ***
  DE2NON:ratioprot -0.04783    0.26268  -0.182    0.856    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4381 on 8545 degrees of freedom
Multiple R-squared:  0.6047,	Adjusted R-squared:  0.6045
###################################

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cdslength_ratio_rand_youngold.pdf", width=8, height=8)
ggplot(diff_prot_length_rand, aes(x=ratiocds, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(0,2) +
  labs(x='Ratio of coding sequence length', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m2 <- lm(abs ~ DE2/ratiocds-1, data = diff_prot_length_rand)
summary(y_m2)
#########################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE           -1.14920    0.26002  -4.420    1e-05 ***
  DE2NON           0.16041    0.23702   0.677    0.499    
DE2DE:ratiocds   3.03097    0.25894  11.705   <2e-16 ***
  DE2NON:ratiocds  0.05361    0.23683   0.226    0.821    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4381 on 8545 degrees of freedom
Multiple R-squared:  0.6047,	Adjusted R-squared:  0.6045 
F-statistic:  3268 on 4 and 8545 DF,  p-value: < 2.2e-16
###########################

#load data on 12feb2019
diff_prot_length <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/diff_protein_length.txt', header = T)
str(diff_prot_length)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_diff_protenlength_ratio_youngold.pdf", width=8, height=8)
ggplot(data=diff_prot_length, aes(x=genomcom,y=prop,fill=factor(DE))) +
  scale_fill_manual(values = c("firebrick3","grey"), labels=c("DE","Non-DE"), name="Bias") + 
  ylim(0,0.8) +
  geom_bar(position="dodge",stat="identity",width=0.6) + 
  geom_text(aes(label=number),position=position_dodge(width=0.6), hjust=1.1) +
  coord_flip() +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Proportion of genes', y='Proportion of genes')
dev.off()

#load propotion data, modified codes on 04.Feb.2019.
stopcodon_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/Mvsl_a1a2_exp_gencompt_protlength_fi.txt', header = T)
str(stopcodon_ratio)

stopcodon_ratio_rmcentro <- subset(stopcodon_ratio,stopcodon_ratio$youngold != "Centro")
str(stopcodon_ratio_rmcentro)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protenlength_ratio_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=ratioprot, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,1.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of protein length (A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cdslength_ratio_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=ratiocds, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.5,1.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length (A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_protlength_corre_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=ratioprot, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-8,12) + xlim(0.5,2) +
  labs(x='Ratio of protein length (A1/A2)', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cdslength_corre_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=ratiocds, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-8,12) + xlim(0.5,2) +
  labs(x='Ratio of coding sequence length (A1/A2)', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <-lm(logFC.A1.A2~ratioprot*DE, data=stopcodon_ratio_rmcentro)
summary(y)
y1 <-lm(logFC.A1.A2~ratioprot+DE, data=stopcodon_ratio_rmcentro)
summary(y1)
anova(y,y1)

y2 <-lm(logFC.A1.A2~ratiocds*DE, data=stopcodon_ratio_rmcentro)
summary(y2)
y3 <-lm(logFC.A1.A2~ratiocds+DE, data=stopcodon_ratio_rmcentro)
summary(y3)
anova(y2,y3)


#stats on 4Feb.2019
stopcodon_oldstrata <- subset(stopcodon_ratio_rmcentro, stopcodon_ratio_rmcentro$youngold=="OldStrata")
wilcox.test(stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='Up'],stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='NON'], exact = FALSE) 
#W = 1811, p-value = 0.1016
wilcox.test(stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='Down'],stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='NON'], exact = FALSE) 
#W = 1856, p-value = 0.1494
wilcox.test(stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='Down'],stopcodon_oldstrata$ratioprot[stopcodon_oldstrata$DE=='Up'], exact = FALSE) 
#W = 661, p-value = 0.8872

gene_ratio_prot_notequal <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/04Feb2019/number_gene_ratio.txt', header = T)
str(gene_ratio_prot_notequal)

gene_ratio_prot_notequal_prot <- subset(gene_ratio_prot_notequal, gene_ratio_prot_notequal$type == "proteinanot")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratioprot_not1_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_prot_notequal_prot, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE","Non-DE"), name="Bias") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  ylim(0,0.1) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion of genes in each category')
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiocds_not1_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_prot_notequal_prot, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE","Non-DE"), name="Bias") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  ylim(0,0.1) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion of genes in each category')
dev.off()


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
  ylim(0.8,1.2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencea2_devidprotein3times_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=cdsa2expest, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.8,1.2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequenceratio_youngold.pdf", width=8, height=8)
ggplot(stopcodon_ratio_rmcentro, aes(x=youngold, y=cdsratioa1a2)) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.9,1.1) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length (A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

stopcodon_ratio_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/19jan2019/Mvsl_a1a2_exp_cds_protein_compart_sep.txt', header = T)
str(stopcodon_ratio_sep)

stopcodon_ratio_sep1 <- subset(stopcodon_ratio_sep, stopcodon_ratio_sep$youngold != "Centro")
str(stopcodon_ratio_sep1)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencea1a2_devidprotein3times_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_sep1, aes(x=youngold, y=cdsa2expest, fill=interaction(haploid,DE))) +
         scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
         geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
         ylim(0.9,1.1) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(cdsa2expest ~ haploid*youngold, data = stopcodon_ratio_sep1)
summary(y)
#########
Coefficients:
  Estimate Std. Error   t value Pr(>|t|)    
(Intercept)                    1.004e+00  4.465e-05 22484.262   <2e-16 ***
  haploidA2                      3.370e-07  6.314e-05     0.005    0.996    
youngoldbPAR                   5.578e-04  3.754e-04     1.486    0.137    
youngoldColorStrata           -5.152e-05  6.676e-04    -0.077    0.938    
youngoldOldStrata             -2.794e-04  2.934e-04    -0.952    0.341    
haploidA2:youngoldbPAR        -2.311e-07  5.308e-04     0.000    1.000    
haploidA2:youngoldColorStrata -3.447e-05  9.442e-04    -0.037    0.971    
haploidA2:youngoldOldStrata    4.505e-05  4.149e-04     0.109    0.914    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.003997 on 16700 degrees of freedom
Multiple R-squared:  0.000364,	Adjusted R-squared:  -5.502e-05 
F-statistic: 0.8687 on 7 and 16700 DF,  p-value: 0.5304
##########

###proportion of genes with length difference, loading on Jan.20.2019.
gene_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/19jan2019/number_gene_ratio.txt', header = T)
str(gene_ratio)

gene_ratio_cds <- subset(gene_ratio, gene_ratio$category =="cds")
gene_ratio_prot <- subset(gene_ratio, gene_ratio$category =="prot")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiogene_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_cds, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","dodgerblue4","light grey","dark grey"), labels=c("DE & ratio!=1","DE & ratio=1","Non-DE & ratio!=1","Non-DE & ratio=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiogeneprot_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_prot, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","dodgerblue4","light grey","dark grey"), labels=c("DE & ratio!=1","DE & ratio=1","Non-DE & ratio!=1","Non-DE & ratio=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

gene_ratio_cds_notequal <- subset(gene_ratio_cds, gene_ratio_cds$type =="cdsnot")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiocsd_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_cds_notequal, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE & ratio!=1","Non-DE & ratio!=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

gene_ratio_prot_notequal <- subset(gene_ratio_prot, gene_ratio_prot$type =="proteinanot")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratioprot_youngold1.pdf", width=8, height=8)
ggplot(gene_ratio_prot_notequal, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE & ratio!=1","Non-DE & ratio!=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
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