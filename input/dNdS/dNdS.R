#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
dNdS <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_TEinsert_2k10kupdownstream_dNdS.txt', header = T)
str(dNdS)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dS_3genomiccompartments.pdf", width=8, height=8)
ggplot(dNdS, aes(x=logFC.A1.A2, y=dS, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(width = .5)) +
  facet_grid(~genomiccomp) +
  ylim(0,0.1) +  
  xlim(-4,4) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dS_8genomiccompartments.pdf", width=8, height=8)
ggplot(dNdS, aes(x=logFC.A1.A2, y=dS, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(width = .5)) +
  facet_grid(~chr) +
  ylim(0,0.1) +  
  xlim(-4,4) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###some stats

y <- lm(logFC.A1.A2~dS*genomiccomp-1, data= dNdS)
summary(y)

####
lm(formula = logFC.A1.A2 ~ dS * genomiccomp - 1, data = dNdS)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.6649  -0.1326   0.0035   0.1382   5.3486 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dS                   -1.747472   1.725681  -1.013  0.31128    
genomiccompaAutosome  0.024075   0.005013   4.802  1.6e-06 ***
  genomiccompbPAR       0.040399   0.043377   0.931  0.35171    
genomiccompNRR        0.116573   0.042438   2.747  0.00603 ** 
  dS:genomiccompbPAR    3.877068  43.781446   0.089  0.92944    
dS:genomiccompNRR     2.425015   2.055512   1.180  0.23814 
####

y1 <- lm(logFC.A1.A2~dS*chr-1, data= dNdS)
summary(y1)
####
lm(formula = logFC.A1.A2 ~ dS * chr - 1, data = dNdS)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.6562  -0.1326   0.0038   0.1385   4.9374 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dS             -1.747472   1.709805  -1.022  0.30681    
chraAutosome    0.024075   0.004967   4.847 1.28e-06 ***
  chrbPAR         0.040399   0.042978   0.940  0.34725    
chrcGreen      -0.217167   0.440557  -0.493  0.62208    
chrdRed         0.062902   0.093284   0.674  0.50014    
chreOrange      3.557926   0.586413   6.067 1.38e-09 ***
  chrfaBlack      0.129097   0.053869   2.396  0.01658 *  
  chrfBlue       -0.209214   0.171284  -1.221  0.22197    
chrgPurple      0.072328   0.195192   0.371  0.71099    
dS:chrbPAR      3.877068  43.378664   0.089  0.92878    
dS:chrcGreen   51.974721  82.969514   0.626  0.53106    
dS:chrdRed     16.500924   7.925997   2.082  0.03740 *  
  dS:chreOrange -35.139061  12.340271  -2.848  0.00442 ** 
  dS:chrfaBlack   1.909564   2.223309   0.859  0.39044    
dS:chrfBlue    -2.555424   4.292616  -0.595  0.55166    
dS:chrgPurple   4.050294   3.554415   1.140  0.25453 
####

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dN_8genomiccompartments.pdf", width=8, height=8)
ggplot(dNdS, aes(x=logFC.A1.A2, y=dN, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(width = .5)) +
  facet_grid(~chr) +
  ylim(0,0.1) +  
  xlim(-4,4) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dN_3genomiccompartments.pdf", width=8, height=8)
ggplot(dNdS, aes(x=logFC.A1.A2, y=dN, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(width = .5)) +
  facet_grid(~genomiccomp) +
  ylim(0,0.05) +  
  xlim(-4,4) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

####
y2 <- lm(logFC.A1.A2~dN*chr-1, data= dNdS)
summary(y2)

lm(formula = logFC.A1.A2 ~ dN * chr - 1, data = dNdS)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.6888  -0.1328   0.0043   0.1396   4.9373 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dN            -3.244e+00  2.588e+00  -1.253   0.2102    
chraAutosome   2.411e-02  4.914e-03   4.906 9.52e-07 ***
  chrbPAR        1.755e-02  4.289e-02   0.409   0.6823    
chrcGreen      1.449e-01  3.689e-01   0.393   0.6944    
chrdRed        6.963e-02  9.017e-02   0.772   0.4401    
chreOrange     5.689e+00  5.650e-01  10.069  < 2e-16 ***
  chrfaBlack     1.022e-01  5.593e-02   1.827   0.0678 .  
chrfBlue      -1.504e+00  1.593e-01  -9.441  < 2e-16 ***
  chrgPurple     6.154e-02  1.703e-01   0.361   0.7179    
dN:chrbPAR     1.573e+02  5.158e+01   3.050   0.0023 ** 
  dN:chrcGreen  -1.480e+02  3.435e+02  -0.431   0.6667    
dN:chrdRed     2.445e+01  1.114e+01   2.195   0.0282 *  
  dN:chreOrange -1.629e+02  2.340e+01  -6.962 3.70e-12 ***
  dN:chrfaBlack  5.527e+00  4.075e+00   1.356   0.1751    
dN:chrfBlue    7.669e+01  8.511e+00   9.011  < 2e-16 ***
  dN:chrgPurple  9.167e+00  6.346e+00   1.444   0.1487  


####

y3 <- lm(logFC.A1.A2~dN*genomiccomp-1, data= dNdS)
summary(y3)
lm(formula = logFC.A1.A2 ~ dN * genomiccomp - 1, data = dNdS)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.8124  -0.1326   0.0041   0.1391   5.3182 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dN                    -3.243652   2.635132  -1.231 0.218398    
genomiccompaAutosome   0.024110   0.005003   4.819 1.47e-06 ***
  genomiccompbPAR        0.017554   0.043663   0.402 0.687662    
genomiccompNRR         0.006655   0.042954   0.155 0.876871    
dN:genomiccompbPAR   157.298380  52.506349   2.996 0.002748 ** 
  dN:genomiccompNRR     13.080650   3.560540   3.674 0.000241 ***
####
  
DE_dNdS_sub <-subset(dNdS, dNdS$genomiccomp == "NRR")
str(DE_dNdS_sub)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dn_onlyNRR.pdf", width=8, height=8)
ggplot(DE_dNdS_sub, aes(x=logFC.A1.A2, y=dN, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA) +
  ylim(0,0.04) +  
  xlim(-4,3) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

wilcox.test(DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='Up'],DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='NON'],exact = FALSE) 
#W = 2352, p-value = 5.645e-05
wilcox.test(DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='Down'],DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='NON'],exact = FALSE) 
#W = 1157.5, p-value = 0.05077
wilcox.test(DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='Down'],DE_dNdS_sub$dN[DE_dNdS_sub$DE_status=='Up'],exact = FALSE) 
#W = 92, p-value = 0.06048

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_ds_onlyNRR.pdf", width=8, height=8)
ggplot(DE_dNdS_sub, aes(x=logFC.A1.A2, y=dS, color=DE_status)) + 
  geom_boxplot(outlier.shape=NA) +
  ylim(0,0.1) +  
  xlim(-4,3) +
  labs(y='dS') + 
  labs(x='LogFC(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

wilcox.test(DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='Up'],DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='NON'],exact = FALSE) 
#W = 2181, p-value = 0.001478
wilcox.test(DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='Down'],DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='NON'],exact = FALSE) 
#W = 1279, p-value = 0.005356
wilcox.test(DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='Down'],DE_dNdS_sub$dS[DE_dNdS_sub$DE_status=='Up'],exact = FALSE) 
#W = 166, p-value = 0.5982


#### between species dNdS.
dNdS_2species <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_a1a2_DEnonDE_betweenspecies_dNdS.txt', header = T)
str(dNdS_2species)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dndifferce_2species.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dndiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("A2 bias", "not bias","A1 bias"), name="DE expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,4) +
  labs(y='dN between Mvls and Mvld difference (A1-A2)') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dsdifferce_2species.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("A2 bias", "not bias","A1 bias"), name="DE expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
 ylim(-2,4) +
  labs(y='dS between Mvls and Mvld difference (A1-A2)') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dndsdifferce_2species.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dndsdiff,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("A2 bias", "not bias","A1 bias"), name="DE expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +
  labs(y='dN/dS between Mvls and Mvld difference (A1-A2)') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
dNdS_2species_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_a1a2_DEnonDE_betweenspecies_dNdS_sep.txt', header = T)
str(dNdS_2species_sep)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_a1a2sep_2species.pdf", width=8, height=8)
ggplot(dNdS_2species_sep, aes(x=chrom, y=dN,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","grey","dark grey","dodgerblue2","dodgerblue4")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,3.5) +
  labs(y='dN between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


dNdS_a2 <- subset(dNdS_2species_sep, dNdS_2species_sep$haploid == "A2")
str(dNdS_a2)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dndsdifferce_2species.pdf", width=8, height=8)
ggplot(dNdS_a2, aes(x=chrom, y=dN,fill=DE)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.02) +
  labs(y='dN between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dS_a1a2sep_2species.pdf", width=8, height=8)
ggplot(dNdS_2species_sep, aes(x=chrom, y=dS,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","grey","dark grey","dodgerblue2","dodgerblue4")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,5) +
  labs(y='dS between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dNdS_a1a2sep_2species.pdf", width=8, height=8)
ggplot(dNdS_2species_sep, aes(x=chrom, y=dNdS,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","pink","grey","white","dodgerblue2","light blue")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,2) +
  labs(y='dN/dS between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###filtering out large dNdS and dS data
dNdS_filter_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_a1a2_DEnonDE_betweenspecies_dNdS_sep copy.txt', header = T)
str(dNdS_filter_sep)

ggplot(dNdS_filter_sep, aes(x=chrom, y=dS,fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","grey","dark grey","dodgerblue2","dodgerblue4")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,3) +
  labs(y='dS between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

