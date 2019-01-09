#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load the corresponding data files, between a1 and a2 homologs within Mvsl species, modified at jan.08.2019
dNdS <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_TEinsert_2k10kupdownstream_dNdS.txt', header = T)
str(dNdS)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_dNdS_3genomiccompartments_v2.pdf", width=8, height=8)
ggplot(dNdS, aes(x=genomiccomp, y=dN, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) +
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
  scale_fill_manual(values = c("firebrick2", "light grey","dodgerblue2")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +
    scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
    labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###some stats

y <- lm(log2exp~dsdiff*chrom-1, data= dNdS)
summary(y)

####
lm(formula = log2exp ~ dsdiff * chrom - 1, data = dNdS)

Residuals:
  Min      1Q  Median      3Q     Max 
-6.0882 -0.0219  0.0163  0.0496  6.4904 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dsdiff            0.016162   0.054505   0.297 0.766841    
chromaAutosome   -0.013779   0.004094  -3.366 0.000768 ***
  chrombPAR         0.014794   0.036146   0.409 0.682349    
chromNRR          0.076732   0.025226   3.042 0.002362 ** 
  dsdiff:chrombPAR -1.173887   0.291709  -4.024 5.79e-05 ***
  dsdiff:chromNRR   0.139897   0.169498   0.825 0.409201    
####

y1 <- lm(log2exp~dndiff*chrom-1, data= dNdS)
summary(y1)
####
lm(formula = log2exp ~ dndiff * chrom - 1, data = dNdS)

Residuals:
  Min      1Q  Median      3Q     Max 
-6.0882 -0.0218  0.0163  0.0496  6.4904 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dndiff            0.054073   0.127718   0.423 0.672037    
chromaAutosome   -0.013782   0.004094  -3.366 0.000766 ***
  chrombPAR         0.014921   0.036149   0.413 0.679787    
chromNRR          0.075833   0.025220   3.007 0.002650 ** 
  dndiff:chrombPAR -2.323792   0.576218  -4.033 5.58e-05 ***
  dndiff:chromNRR   0.363448   0.421854   0.862 0.388970  

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
y2 <- lm(log2exp~dndsdiff*chrom-1, data= dNdS)
summary(y2)

###
  Call:
  lm(formula = log2exp ~ dndsdiff * chrom - 1, data = dNdS)

Residuals:
  Min      1Q  Median      3Q     Max 
-6.0881 -0.0214  0.0166  0.0497  6.4905 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
dndsdiff            6.125e-11  9.602e-11   0.638 0.523607    
chromaAutosome     -1.388e-02  4.102e-03  -3.383 0.000721 ***
  chrombPAR          -1.544e-02  3.716e-02  -0.416 0.677773    
chromNRR            6.100e-02  2.666e-02   2.288 0.022144 *  
  dndsdiff:chrombPAR  9.258e-03  6.641e-03   1.394 0.163386    
dndsdiff:chromNRR   8.540e-11  1.269e-10   0.673 0.500979 

####
####
###dn stats in NRR
wilcox.test(dNdS$dndiff[dNdS$DE=='Up'],dNdS$dndiff[dNdS$DE=='NON'],exact = FALSE) 
#W = 270080, p-value = 9.267e-09
wilcox.test(dNdS$dndiff[dNdS$DE=='Down'],dNdS$dndiff[dNdS$DE=='NON'],exact = FALSE)  
#W = 276300, p-value = 0.2152
wilcox.test(dNdS$dndiff[dNdS$DE=='Up'],dNdS$dndiff[dNdS$DE=='Down'],exact = FALSE) 
#W = 4046.5, p-value = 0.08721


###dS stats in NRR
wilcox.test(dNdS$dsdiff[dNdS$DE=='Up'],dNdS$dsdiff[dNdS$DE=='NON'],exact = FALSE) 
#W = 272720, p-value = 6.103e-05
wilcox.test(dNdS$dsdiff[dNdS$DE=='Down'],dNdS$dsdiff[dNdS$DE=='NON'],exact = FALSE)  
#W = 279030, p-value = 0.2809
wilcox.test(dNdS$dsdiff[dNdS$DE=='Up'],dNdS$dsdiff[dNdS$DE=='Down'],exact = FALSE) 
#W = 4053.5, p-value = 0.1041


###dn/dS stats in NRR
wilcox.test(dNdS$dndsdiff[dNdS$DE=='Up'],dNdS$dndsdiff[dNdS$DE=='NON'],exact = FALSE) 
#W = 291950, p-value = 0.0002242
wilcox.test(dNdS$dndsdiff[dNdS$DE=='Down'],dNdS$dndsdiff[dNdS$DE=='NON'],exact = FALSE)  
#W = 326470, p-value = 0.0003383
wilcox.test(dNdS$dndsdiff[dNdS$DE=='Up'],dNdS$dndsdiff[dNdS$DE=='Down'],exact = FALSE) 
#W = 3834.5, p-value = 0.5456


#### between 2 sub-species dNdS.
dNdS_2species <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/dNdS/Mvsl_DEnonDE_2subspecies_dnds_sep_fi.txt', header = T)
str(dNdS_2species)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dn_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dN, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.015) +
  labs(y='dN between MvSl-1064 and MvSl-1318') +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) +
  labs(x='Genomic compartment') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_dS_2subspecies_3compartm.pdf", width=8, height=8)
ggplot(dNdS_2species, aes(x=chrom, y=dS, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) +
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
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) + 
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
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
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
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
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
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
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
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) + 
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
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) + 
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
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) + 
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
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2")) + 
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
