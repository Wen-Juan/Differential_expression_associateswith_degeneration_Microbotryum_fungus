#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#loading data on 01Feb.2019.
##below are genes with single copy homologs.
TE_number_touse_singlecopy_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/TE_number_touse_01feb19.txt', header = T)
str(TE_number_touse_singlecopy_homolog)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_prop_all_intervals_01feb19.pdf", width=8, height=8)
ggplot(TE_number_touse_singlecopy_homolog, aes(x=location,y=Te_prop, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/geneswithte_prop_all_intervals_01feb19.pdf", width=8, height=8)
ggplot(TE_number_touse_singlecopy_homolog, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of genes with TE insertions', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### loading dataset of individual gene with TE insertions, 03Feb.2019
TE_singlegene <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/Mvsl_a1a2_te_exp_compart.txt', header = T)
str(TE_singlegene)

TE_singlegene_rmcentro <- subset(TE_singlegene, TE_singlegene$youngold!="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genediff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=genediff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Gene difference in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_geneup5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down5kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown10kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down10kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown15kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down15kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown20kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down20kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up10kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 5-10k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up15kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 10-15k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up20kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 15-20k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_gene0k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=genediff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for homologs introns (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down5k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down5kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 0-5k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down5k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down5kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 0-5k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up10k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk10diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 5-10k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down10k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down10kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 5-10k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down15k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down15kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 10-105k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down20k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down20kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 15-20k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up15k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk15diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 10-15k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up20k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk20diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 15-20k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###stats.
y <- lm(logFC.A1.A2 ~ upk5diff*DE, data = TE_singlegene_rmcentro)
summary(y)
y1 <- lm(logFC.A1.A2 ~ upk5diff+DE, data = TE_singlegene_rmcentro)
anova(y,y1)
###
Model 1: logFC.A1.A2 ~ upk5diff * DE
Model 2: logFC.A1.A2 ~ upk5diff + DE
Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
1   8542 2022.3                                
2   8544 2025.0 -2   -2.7503 5.8085 0.003014 **
###
 
y2 <- lm(logFC.A1.A2 ~ upk10diff*DE, data = TE_singlegene_rmcentro)
summary(y2)
y3 <- lm(logFC.A1.A2 ~ upk10diff+DE, data = TE_singlegene_rmcentro)
anova(y2,y3) 
###

y6 <- lm(logFC.A1.A2 ~ down5kdiff*DE, data = TE_singlegene_rmcentro)
summary(y6)
y7 <- lm(logFC.A1.A2 ~ down5kdiff+DE, data = TE_singlegene_rmcentro)
anova(y6,y7) 

TE_singlegene_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/Mvsl_a1a2_te_exp_compart_sep.txt', header = T)
str(TE_singlegene_sep)

TE_singlegene_sep_rmcentro <- subset(TE_singlegene_sep, TE_singlegene_sep$youngold != "Centro")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_geneup5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_sep_rmcentro, aes(x=youngold, y=genea10k, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(logFC.A1.A2 ~ genediff*DE, data = TE_singlegene_rmcentro)
summary(y)
y1 <- lm(logFC.A1.A2 ~ genediff+DE, data = TE_singlegene_rmcentro)
anova(y,y1)
###

##below are all genes
TE_number_touse <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/TE_number_touse.txt', header = T)
str(TE_number_touse)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_prop_all_intervals.pdf", width=8, height=8)
ggplot(TE_number_touse, aes(x=location,y=Te_prop, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.2) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/genes_withTEinsertions_intervals.pdf", width=8, height=8)
ggplot(TE_number_touse, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.2) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#loading data on 24Jan.2019, check distribution of coding sequence length
a1cds_length <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/A1_cds_length.txt', header = F)
str(a1cds_length)
a2cds_length <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/A2_cds_length.txt', header = F)
str(a2cds_length)
length_diff <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/homolog_length_variation.txt', header = F)
str(length_diff)

ggplot(data=a1cds_length, aes(a1cds_length$V2)) + 
  xlim(0,7500) +
  ylim(0,2000) +
  geom_histogram(binwidth=200)

ggplot(data=a2cds_length, aes(a2cds_length$V2)) + 
  xlim(0,7500) +
  geom_histogram(binwidth=200)

ggplot(length_diff, aes(length_diff$V8)) + 
  xlim(0,1000) +
  ylim(0,40) +
  geom_histogram(binwidth=25)

#load the corresponding data files. modify this on Jan.17.2019.
#summary of TE numbers
TE_number <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/TE_nr_17012019.txt', header = T)
str(TE_number)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_number_intervals_5k.pdf", width=8, height=8)
ggplot(TE_number, aes(x=location,y=TE_interval, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,11000) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Number of detected TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Prop_geneswith_TEinsertion_intervals_5k.pdf", width=8, height=8)
ggplot(TE_number, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.32) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of total genes with TE insertions', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#load the corresponding data files, DE genes and TE intervals. modify this on Jan.17.2019.
TE_DEexp <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/Mvsl_TE_exp_allintervals.txt', header = T)
str(TE_DEexp)

TE_DEexp1 <- subset(TE_DEexp, TE_DEexp$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm (logFC.A1.A2~upk20diff*DE-1, data=TE_DEexp1)
summary(y)

###
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
upk5diff       -0.159073   0.041431  -3.840 0.000124 ***
  DEDown         -1.851788   0.031855 -58.132  < 2e-16 ***
  DENON          -0.004183   0.004791  -0.873 0.382622    
DEUp            1.819751   0.036815  49.430  < 2e-16 ***
  upk5diff:DENON  0.192189   0.042263   4.547 5.51e-06 ***
  upk5diff:DEUp   0.237904   0.051473   4.622 3.86e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4227 on 8348 degrees of freedom
Multiple R-squared:  0.4217,	Adjusted R-squared:  0.4213 
F-statistic:  1015 on 6 and 8348 DF,  p-value: < 2.2e-16
###

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5-10kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5-10kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_10-15kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_10-15kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_15-20kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_15-20kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###only working on data on MAT chromosome.
TE_DEexp_mat <- subset(TE_DEexp, TE_DEexp$youngold != "Auto")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_0-5kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk5diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() +geom_smooth(method=lm) +
  ylim(-12,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 5k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_0-5kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk5diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-6,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 5k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_5-10kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk10diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 5-10k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_5-10kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk10diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 5-10k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_10-15kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk15diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 10-15k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_10-15kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk15diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 10-15k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_15-20kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk20diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 15-20k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_15-20kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk20diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 15-20k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load TE and exp. data separately, Jan.18.2018
UpTE_DEexp_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/UpstreamTE_Exp_intervals_sep.txt', header = T)
str(UpTE_DEexp_sep)

UpTE_DEexp_sep1 <- subset(UpTE_DEexp_sep, UpTE_DEexp_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_upstream_exp.pdf", width=8, height=8)
ggplot(UpTE_DEexp_sep1, aes(x=youngold, y=up5k, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +
  facet_grid(~Up) +
  scale_x_discrete(labels=c("Auto", "PAR", "Young","Old")) + 
  labs(x='Genomic compartment', y='Average TE insertion sites') +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

DownTE_DEexp_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/DownstreamTE_Exp_intervals_sep.txt', header = T)
str(DownTE_DEexp_sep)
DownTE_DEexp_sep1 <- subset(DownTE_DEexp_sep, DownTE_DEexp_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_downstream_exp.pdf", width=8, height=8)
ggplot(DownTE_DEexp_sep1, aes(x=youngold, y=downak5, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +
  facet_grid(~Down) +
  scale_x_discrete(labels=c("Auto", "PAR", "Young","Old")) + 
  labs(x='Genomic compartment', y='TE insertion at 5k-upstream (A1-A2)') +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#load the corresponding data files. modify this on Jan.09.2019.
#summary of TE numbers
TE_number <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/TE_number.txt', header = T)
str(TE_number)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_number_intervals.pdf", width=8, height=8)
ggplot(TE_number, aes(x=range, y=TE, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,4000) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Number of detected TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#prop. gene with TE insertions
Prop_gene <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/gene_TEinsertion_prop.txt', header = T)
str(Prop_gene)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Gene_proportion_intervals.pdf", width=8, height=8)
ggplot(Prop_gene, aes(x=range, y=proportion, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#de genes on MAT 
exp_DE_MAT_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_mat_cor.txt', header = T)
str(exp_DE_MAT_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_diffMAT_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_MAT_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.25,0.4) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

dev.off()


#DE autosome
exp_DE_auto_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_auto_cor.txt', header = T)
str(exp_DE_auto_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_autosomes_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_auto_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.2,0.2) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#################################not consider this at the moment################################################
### richmedium
DE_richmedium <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_richmedium.txt', header = T)
str(DE_richmedium)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDEinterval_overlapwithTE_propgenes_richmeidum_cor.pdf", width=8, height=8)
ggplot(DE_richmedium, aes(x=interval, y=diffpropcor, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.1,0.1) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
#############################################################################################

### 1,20k up and down stream TE insertion sites among genomic compartments.
DE_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_2k10kupdownstream.txt', header = T)
str(DE_homolog)

################### I think this figure is not correct, both catogory makes the boxplot look difference among DE and Non-DE genes.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog, aes(y=logFC.A1.A2, x=k2updiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~chr) +
  ylim(-4,4) +  
  xlim(-2,2) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
ddev.off()
############################################################################################################
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff2kdownstream_logfc_3genomicparts.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-0.5,0.5) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff2kupstream_logfc_3genomicparts,.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-0.5,0.5) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-2k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff8kupstream_logfc_3genomicparts", width=8, height=8)
  ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
    labs(x='Genomic compartment') + 
    labs(y='Difference of TE number [2-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
  
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kdwonstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervaldowndiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
DE_homolog2 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_mod8.txt', header = T)
str(DE_homolog2)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20updiffa1,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20downdiffa2,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

