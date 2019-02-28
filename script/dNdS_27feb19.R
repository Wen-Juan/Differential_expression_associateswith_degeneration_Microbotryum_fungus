#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
install.packages("gridExtra")
library(gridExtra)
require(cowplot)

#load the corresponding data files, between a1 and a2 homologs within Mvsl species, homolog length was trimmed to equal and removing TE genes. modified at jan.28.2019
dNdS_new <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/28jan2019/Mvsl_dnds_exp_gencomp.txt', header = T)
str(dNdS_new)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dN_DEnonDe_twoparts28012019.pdf", width=8, height=8)
p_dN <- ggplot(dNdS_new, aes(x=youngold, y=dn, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ds_DEnonDe_twoparts28012019.pdf", width=8, height=8)
p_dS <- ggplot(dNdS_new, aes(x=youngold, y=ds, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

#combine two figures in one
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_ds_combine_new.pdf", width=8, height=8)
plot_grid(p_dN, p_dS, labels = c('A', 'B'))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dnds_DEnonDe_twoparts28012019.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=youngold, y=dnds, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,7) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


##stats on 28jan2019
dNdS_youngstrata <- subset(dNdS_new, dNdS_new$youngold=="ColorStrata")
dNdS_oldstrata <- subset(dNdS_new, dNdS_new$youngold=="OldStrata")
  
wilcox.test(dNdS_youngstrata$dn[dNdS_youngstrata$DE2=='DE'],dNdS_youngstrata$dn[dNdS_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 74, p-value = 0.0331
wilcox.test(dNdS_youngstrata$ds[dNdS_youngstrata$DE2=='DE'],dNdS_youngstrata$ds[dNdS_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 75, p-value = 0.02791
wilcox.test(dNdS_youngstrata$dnds[dNdS_youngstrata$DE2=='DE'],dNdS_youngstrata$dnds[dNdS_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 47, p-value = 0.7608

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 2941, p-value = 0.0003843
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$ds[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 2960.5, p-value = 0.000277
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 2276, p-value = 0.5286

#stats, separating A1- and A2-biased genes.
dNdS_youngstrata <- subset(dNdS_new, dNdS_new$youngold=="ColorStrata")
dNdS_oldstrata <- subset(dNdS_new, dNdS_new$youngold=="OldStrata")

wilcox.test(dNdS_youngstrata$ds[dNdS_youngstrata$DE=='Down'],dNdS_youngstrata$ds[dNdS_youngstrata$DE=='NON'],exact = FALSE) 
#W = 75, p-value = 0.02791
wilcox.test(dNdS_youngstrata$dn[dNdS_youngstrata$DE=='Down'],dNdS_youngstrata$dn[dNdS_youngstrata$DE=='NON'],exact = FALSE) 
#W = 74, p-value = 0.0331
wilcox.test(dNdS_youngstrata$dnds[dNdS_youngstrata$DE=='Down'],dNdS_youngstrata$dnds[dNdS_youngstrata$DE=='NON'],exact = FALSE) 
#W = 47, p-value = 0.7608

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1143.5, p-value = 0.0163
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1020.5, p-value = 0.1466

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1797.5, p-value = 0.003051
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1940, p-value = 0.0001384

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 211, p-value = 0.779
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 142, p-value = 0.1244

wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1077, p-value = 0.05896
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1199, p-value = 0.5491
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 287, p-value = 0.02079

###scatter point figures, on 29.jan.2019
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_exp_correlation_DEpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=dn, y=abs,color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_exp_correlation_notpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=dn, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ds_exp_correlation_pool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=ds, y=abs,color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ds_exp_correlation_DEnotpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=ds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dnds_exp_correlation_DEpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=dnds, y=abs,color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("Biased","Not-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(0,13) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dnds_exp_correlation_DEnotpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=dnds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(0,13) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(abs ~ dn*DE2, data = dNdS_new)
summary(y)
y1 <- lm(abs ~ dn+DE2, data = dNdS_new)
summary(y1)
anova(y,y1)
###
Model 1: abs ~ dn * DE2
Model 2: abs ~ dn + DE2
Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
1   6110 347.42                                  
2   6111 358.30 -1   -10.887 191.47 < 2.2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###

y <- lm(abs ~ dn*DE, data = dNdS_new)
summary(y)
y1 <- lm(abs ~ dn+DE, data = dNdS_new)
summary(y1)
anova(y,y1)
###
Analysis of Variance Table

Model 1: abs ~ dn * DE
Model 2: abs ~ dn + DE
Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
1   6108 343.47                                  
2   6110 355.42 -2   -11.953 106.28 < 2.2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###
###

y2 <- lm(abs ~ ds*DE2, data = dNdS_new)
summary(y2)
y3 <- lm(abs ~ ds+DE2, data = dNdS_new)
summary(y3)
anova(y2,y3)
###
Analysis of Variance Table

Model 1: abs ~ ds * DE2
Model 2: abs ~ ds + DE2
Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
1   6110 344.79                                  
2   6111 360.18 -1   -15.392 272.76 < 2.2e-16 ***
###

y4 <- lm(abs ~ dnds*DE2, data = dNdS_new)
summary(y4)
y5 <- lm(abs ~ dnds+DE2, data = dNdS_new)
summary(y5)
anova(y4,y5)

###
Analysis of Variance Table

Model 1: abs ~ dnds * DE2
Model 2: abs ~ dnds + DE2
Res.Df    RSS Df Sum of Sq      F Pr(>F)
1   6110 373.07                           
2   6111 373.12 -1 -0.047997 0.7861 0.3753
###

##load data of dN/dS between Mvsl and Mvsv, on 31.Jan.2019.
dNdS_MvslMvsv <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/28jan2019/Mvsl_Mvsv/Mvsl_Mvsv_dnds_exp_all_fi.txt', header = T)
str(dNdS_MvslMvsv)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_dNdsdiffe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dndsdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.3,0.3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN/dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_dNdiffe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dndiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.03,0.03) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_dSdiffe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dsdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.04,0.04) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

dNdS_MvslMvsv_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/28jan2019/Mvsl_Mvsv/Mvsl_Mvsv_dnds_exp_all_fi_sep.txt', header = T)
str(dNdS_MvslMvsv_sep)
dNdS_MvslMvsv_sep_oldstrata <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold=="OldStrata")
x <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Lowmutations")
y <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Vhighmutations")
z <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Neutral")

mean(x$dna2) #0.03965833
mean(x$dsa2) #0.1263271
mean(x$dndsa2) # 0.3424813

mean(z$dna2) #0.03380762
mean(z$dsa2) #0.1370848
mean(z$dndsa2) #0.2839695

mean(y$dna2) #0.0394625
mean(y$dsa2) # 0.1274812
mean(y$dndsa2) # 0.336725

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_dN_mutations_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dna2, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("Low mutations","Equal mutations","High mutations"), name="Expectations") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.45) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_ds_mutations_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dsa2, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("Low mutations","Equal mutations","High mutations"), name="Expectations") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.45) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsv_dNds_mutations_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dndsa2, fill=DE2)) + 
 scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("Low mutations","Equal mutations","High mutations"), name="Expectations") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###stats 31.jan.2019
dNdS_MvslMvsv_sep1 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "OldStrata")

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1165, p-value = 0.927
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1167.5, p-value = 0.9125
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1155.5, p-value = 0.9825

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 4328, p-value = 0.1272
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5302.5, p-value = 0.5743
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5748, p-value = 0.1293

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 4227, p-value = 0.08152
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5240.5, p-value = 0.6681
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5712, p-value = 0.15

dNdS_MvslMvsv_sep2 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "ColorStrata")
wilcox.test(dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],exact = FALSE) 
#W = 5, p-value = 1
wilcox.test(dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],exact = FALSE) 
#W = 4, p-value = 1
wilcox.test(dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],exact = FALSE) 
#W = 5, p-value = 1

wilcox.test(dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 70, p-value = 0.5821
wilcox.test(dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 121, p-value = 0.2638
wilcox.test(dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 132, p-value = 0.1377

wilcox.test(dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],dNdS_MvslMvsv_sep2$dsa2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 69, p-value = 0.5594
wilcox.test(dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],dNdS_MvslMvsv_sep2$dna2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 123, p-value = 0.2364
wilcox.test(dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Lowmutations'],dNdS_MvslMvsv_sep2$dndsa2[dNdS_MvslMvsv_sep2$DE2=='Neutral'],exact = FALSE) 
#W = 133, p-value = 0.1291

dNdS_MvslMvsv_sep3 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "Auto")
wilcox.test(dNdS_MvslMvsv_sep3$dsa2[dNdS_MvslMvsv_sep3$DE2=='Lowmutations'],dNdS_MvslMvsv_sep3$dsa2[dNdS_MvslMvsv_sep3$DE2=='Vhighmutations'],exact = FALSE) 

dNdS_MvslMvsv_sep4 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "bPAR")
wilcox.test(dNdS_MvslMvsv_sep4$dsa2[dNdS_MvslMvsv_sep4$DE2=='Lowmutations'],dNdS_MvslMvsv_sep4$dsa2[dNdS_MvslMvsv_sep4$DE2=='Vhighmutations'],exact = FALSE) 


##laod data of dN/dS between Mvsl and Mvsd, on 30.Jan.2019.
dNdS_MvslMvsd <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/28jan2019/Mvsl_Mvsd/Mvsl_Mvsd_dnds_exp_genomcompt.txt', header = T)
str(dNdS_MvslMvsd)

dNdS_MvslMvsd_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/28jan2019/Mvsl_Mvsd/Mvsl_Mvsd_dnds_exp_genomcompt_sep.txt', header = T)
str(dNdS_MvslMvsd_sep)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dNdiff_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd, aes(x=youngold, y=dndiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.01,0.01) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dsdiff_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd, aes(x=youngold, y=dsdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.025,0.025) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

dNdS_MvslMvsd_oldstrata <- subset(dNdS_MvslMvsd, dNdS_MvslMvsd$youngold == "OldStrata")
str(dNdS_MvslMvsd_oldstrata)

wilcox.test(dNdS_MvslMvsd_oldstrata$dsdiff[dNdS_MvslMvsd_oldstrata$DE=='Up'],dNdS_MvslMvsd_oldstrata$ds[dNdS_MvslMvsd_oldstrata$DE=='NON'],exact = FALSE) 
#V = 78, p-value = 0.9622
wilcox.test(dNdS_MvslMvsd_oldstrata$dsdiff[dNdS_MvslMvsd_oldstrata$DE=='Down'],dNdS_MvslMvsd_oldstrata$ds[dNdS_MvslMvsd_oldstrata$DE=='NON'],exact = FALSE) 
#V = 195.5, p-value = 0.8733
wilcox.test(dNdS_MvslMvsd_oldstrata$dsdiff[dNdS_MvslMvsd_oldstrata$DE=='Down'],dNdS_MvslMvsd_oldstrata$ds[dNdS_MvslMvsd_oldstrata$DE=='Up'],exact = FALSE) 
#V = 195.5, p-value = 0.8733

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dndsdiff_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd, aes(x=youngold, y=dndsdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-1.2,1.2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN/dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_mutationsdn_pool_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd_sep, aes(x=youngold, y=dna1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dodgerblue3","light grey"),labels=c("More mutation","Low mutation","Equal mutation"), name="Mutations bet. homologs") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.03) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_mutationsds_pool_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd_sep, aes(x=youngold, y=dsa1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dodgerblue3","light grey"),labels=c("More mutation","Low mutation","Equal mutation"), name="Mutations bet. homologs") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.04) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_mutationsdnds_pool_youngandold.pdf", width=8, height=8)
ggplot(dNdS_MvslMvsd_sep, aes(x=youngold, y=dndsa1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dodgerblue3","light grey"),labels=c("More mutation","Low mutation","Equal mutation"), name="Mutations bet. homologs") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


ggplot(dNdS_MvslMvsd_oldstrata, aes(x=dsdiff, y=logFC.A1.A2,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load the corresponding data files, between a1 and a2 homologs within Mvsl species, modified at jan.15.2019
dNdS_70perc <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/15jan19/LogCPM_0.05_hwseventyperc_touse_dnds_rmTE.txt', header = T)
str(dNdS_70perc)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dN_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_70perc, aes(x=youngold, y=dn, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dS_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_70perc, aes(x=youngold, y=ds, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dNdS_DEnonDe_youngandold.pdf", width=8, height=8)
ggplot(dNdS_70perc, aes(x=youngold, y=dnds, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dN_exp_correlation.pdf", width=8, height=8)
dNdS_70perc$abs <- abs(dNdS_70perc$logFC.A1.A2)
ggplot(dNdS_70perc, aes(x=dn, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_70perc, aes(x=ds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dndS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_70perc, aes(x=dnds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(abs ~ (dn + ds)*DE, data = dNdS_70perc)
summary(y)
#########
Estimate Std. Error t value Pr(>|t|)    
(Intercept)   1.24448    0.03803  32.722  < 2e-16 ***
  dn          -15.85609    3.18112  -4.984 6.39e-07 ***
  ds           70.93732    6.65037  10.667  < 2e-16 ***
  DENON        -1.06266    0.03820 -27.821  < 2e-16 ***
  DEUp          0.11449    0.04931   2.322 0.020281 *  
  dn:DENON     16.56885    3.32187   4.988 6.28e-07 ***
  dn:DEUp       5.14516    3.91732   1.313 0.189085    
ds:DENON    -69.20246    6.82753 -10.136  < 2e-16 ***
  ds:DEUp     -28.17574    7.99132  -3.526 0.000425 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2719 on 6115 degrees of freedom
Multiple R-squared:  0.4204,	Adjusted R-squared:  0.4197 
F-statistic: 554.5 on 8 and 6115 DF,  p-value: < 2.2e-16
###########


###some stats, 15Jan.2019
dnds_oldstrata <-subset(dNdS_70perc,dNdS_70perc$youngold == "OldStrata")
str(dnds_oldstrata)
wilcox.test(dnds_oldstrata$ds[dnds_oldstrata$DE=='Up'],dnds_oldstrata$ds[dnds_oldstrata$DE=='NON'],exact = FALSE) 
#W = 2101.5, p-value = 2.757e-05
wilcox.test(dnds_oldstrata$ds[dnds_oldstrata$DE=='Down'],dnds_oldstrata$ds[dnds_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1153.5, p-value = 0.008051
wilcox.test(dnds_oldstrata$ds[dnds_oldstrata$DE=='Up'],dnds_oldstrata$ds[dnds_oldstrata$DE=='Down'],exact = FALSE) 
#W = 248, p-value = 0.4287

wilcox.test(dnds_oldstrata$dn[dnds_oldstrata$DE=='Up'],dnds_oldstrata$dn[dnds_oldstrata$DE=='NON'],exact = FALSE) 
#W = 2016, p-value = 0.0002169
wilcox.test(dnds_oldstrata$dn[dnds_oldstrata$DE=='Down'],dnds_oldstrata$dn[dnds_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1208.5, p-value = 0.002046
wilcox.test(dnds_oldstrata$dn[dnds_oldstrata$DE=='Up'],dnds_oldstrata$dn[dnds_oldstrata$DE=='Down'],exact = FALSE) 
#W = 206, p-value = 0.8113

dnds_youngstrata <-subset(dNdS_70perc,dNdS_70perc$youngold == "ColorStrata")
str(dnds_youngstrata)
wilcox.test(dnds_youngstrata$dn[dnds_youngstrata$DE=='Up'],dnds_youngstrata$dn[dnds_youngstrata$DE=='NON'],exact = FALSE) 
#W = 83, p-value = 0.1268
wilcox.test(dnds_youngstrata$ds[dnds_youngstrata$DE=='Up'],dnds_youngstrata$ds[dnds_youngstrata$DE=='NON'],exact = FALSE) 
#W = 94, p-value = 0.03071

#loading MvSl 2sub species dNdS data on Jan.16.2019
dNdS_2sub <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/15jan19/Mvsl_2sub_a1a2_dnds_expdata_genomcomp.txt', header = T)
str(dNdS_2sub)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2sub_dN_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2sub, aes(x=youngold, y=dn, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.015) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2sub_dS_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2sub, aes(x=youngold, y=ds, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.03) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2sub_dNdS_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2sub, aes(x=youngold, y=dnds, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,2.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###loading dN, dS data between MvSl and MvSd Jan.16.2018
dNdS_2species <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/15jan19/Mvsl_Mvsd_a1a2_dnds_expdata.txt', header = T)
str(dNdS_2species)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2species_dN_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=youngold, y=dn, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.02) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2species_dS_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=youngold, y=ds, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.04) +  
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2species_dNdS_exp_youngoldstrata.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=youngold, y=dnds, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,2) +  
  #scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###scatter plots
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a1_dN_exp_correlation.pdf", width=8, height=8)
dNdS_2species$abs <- abs(dNdS_2species$logFC.A1.A2)
dNdS_2species_a1 <- subset(dNdS_2species, dNdS_2species$haploid == "A1")
ggplot(dNdS_2species_a1, aes(x=dn, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,0.08) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

dNdS_2species_a2 <- subset(dNdS_2species, dNdS_2species$haploid == "A2")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a2_dN_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_2species_a2, aes(x=dn, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,0.08) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a1_dS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_2species_a1, aes(x=ds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,0.1) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a2_dS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_2species_a2, aes(x=ds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,0.1) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a1_dNdS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_2species_a1, aes(x=dnds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,2.5) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_2speceis_a2_dNdS_exp_correlation.pdf", width=8, height=8)
ggplot(dNdS_2species_a2, aes(x=dnds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  xlim(0,2.5) +
  geom_point() + geom_smooth(method = lm) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load the corresponding data files, between a1 and a2 homologs within Mvsl species, modified at jan.08.2019
dNdS <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_TEinsert_2k10kupdownstream_dNdS.txt', header = T)
str(dNdS)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dNdS_3genomiccompartments_v2.pdf", width=8, height=8)
ggplot(dNdS, aes(x=genomiccomp, y=dN, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.04) +  
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_ds_3genomiccompartments.pdf", width=8, height=8)
  ggplot(dNdS, aes(x=genomiccomp, y=dS, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2", "light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +
    scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
    labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


####

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dNds_3genomiccompartments.pdf", width=8, height=8)
ggplot(dNdS, aes(x=log2exp, y=dndsdiff, color=DE)) + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~chr) +
  ylim(0,0.005) +  
  xlim(-1,1) +
  labs(y='dN/dS difference (A1-A2)') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

####
###dn stats in NRR
wilcox.test(dNdS$dN[dNdS$DE_status=='Up'],dNdS$dN[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 293270, p-value < 2.2e-16
wilcox.test(dNdS$dN[dNdS$DE_status=='Down'],dNdS$dN[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 209680, p-value < 2.2e-16
wilcox.test(dNdS$dN[dNdS$DE_status=='Down'],dNdS$dN[dNdS$DE_status=='Up'],exact = FALSE)
#W = 1914.5, p-value = 0.1353


###dS stats in NRR
wilcox.test(dNdS$dS[dNdS$DE_status=='Up'],dNdS$dS[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 290370, p-value < 2.2e-16
wilcox.test(dNdS$dS[dNdS$DE_status=='Down'],dNdS$dS[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 210070, p-value < 2.2e-16
wilcox.test(dNdS$dS[dNdS$DE_status=='Down'],dNdS$dS[dNdS$DE_status=='Up'],exact = FALSE)
#W = 2012, p-value = 0.346



###dn stats in NRR
wilcox.test(dNdS$dN[dNdS$DE_status=='Up'],dNdS$dN[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 293270, p-value < 2.2e-16
wilcox.test(dNdS$dN[dNdS$DE_status=='Down'],dNdS$dN[dNdS$DE_status=='NON'],exact = FALSE) 
#W = 209680, p-value < 2.2e-16
wilcox.test(dNdS$dN[dNdS$DE_status=='Down'],dNdS$dN[dNdS$DE_status=='Up'],exact = FALSE)
#W = 1914.5, p-value = 0.1353


#### between 2 sub-species dNdS.
dNdS_2species <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_2subspecies_dnds_sep_fi.txt', header = T)
str(dNdS_2species)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dN, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.015) +
  labs(y='dN between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

ggplot(dNdS_2species, aes(x=dN, y=expmean),color=haploid) +
         #, shape=haploid, color=DE)) +
         geom_point() + geom_smooth(method = lm)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dS_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dS, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.03) +
  labs(y='dS between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dnds_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dnds,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.003) +
  labs(y='dN/dS between MvSl-1064 and MvSl-1318') +
  #scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###difference of dN and dS between 2 sub-species
dNdS_diff <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_2subspecies_dnds_fi_mod.txt', header = T)
str(dNdS_diff)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_diff, aes(x=chrom, y=dndiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.004,0.004) +
  labs(y='dN difference between MvSl-1064 and MvSl-1318') +
  #scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ds_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_diff, aes(x=chrom, y=dsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.004,0.004) +
  labs(y='dS difference between MvSl-1064 and MvSl-1318') +
  #scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dnds_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_diff, aes(x=chrom, y=dndsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.004,0.004) +
  labs(y='dN/dS difference between MvSl-1064 and MvSl-1318') +
  #scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### between Mvsl and Mvsd
dNdS_Mvsl_Mvsd <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_Mvsd/Mvsl_Mvsd_dnds_a1_a2_sep_fi.txt', header = T)
str(dNdS_Mvsl_Mvsd)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dn_a1a2sep_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdS_Mvsl_Mvsd, aes(x=chrom, y=dNa1,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.02) +
  labs(y='dN between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dS_a1a2sep_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdS_Mvsl_Mvsd, aes(x=chrom, y=dSa1,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.04) +
  labs(y='dS between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dNdS_a1a2sep_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdS_Mvsl_Mvsd, aes(x=chrom, y=dNdSa1,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1.2) +
  labs(y='dN/dS between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### dN,dS,dN/dS difference
dNdSdiff_Mvsl_Mvsd <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_Mvsd/Mvsl_Mvsd_dnds_a1a2_match_fi_090119.txt', header = T)
str(dNdSdiff_Mvsl_Mvsd)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dNdiff_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdSdiff_Mvsl_Mvsd, aes(x=chrom, y=dndiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.001,0.001) +
  labs(y='dN difference between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dSdiff_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdSdiff_Mvsl_Mvsd, aes(x=chrom, y=dsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.001,0.001) +
  labs(y='dS difference between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvsd_dndsdiff_2species_3compartm.pdf", width=8, height=8)
ggplot(dNdSdiff_Mvsl_Mvsd, aes(x=chrom, y=dndsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.0005,0.0005) +
  labs(y='dN/dS difference between MvSl and MvSd') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
