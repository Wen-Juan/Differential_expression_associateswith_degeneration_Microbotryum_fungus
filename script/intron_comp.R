#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
intron_total <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/Mvsl_A1A2_singleortholog_genomiccompartment.txt', header = T)
str(intron_total)

intron_diff <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/Mvsl_singleortholog_introndifference.txt', header = T)
str(intron_diff)

#plot intron number
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_introntotallengthovergenelength_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronovergenelength, fill=Haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.75) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron number normalized by coding gene length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#plot intron mean length
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_intronmeanlength_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronmean, fill=Haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,210) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Average intron length per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#plot intron total length
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_introntotallength_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=introntotal, fill=Haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1200) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron total length per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
y1 <- lm(introntotal~Genomicloc, data=intron_total)
y2 <- lm(introntotal~Genomicloc-1, data=intron_total)
anova(y1,y2)
summary(y2)

#Coefficients:
Estimate Std. Error t value Pr(>|t|)    
GenomiclocaAutosome  339.852      2.082 163.202  < 2e-16 ***
  GenomiclocbPAR       304.128     16.904  17.991  < 2e-16 ***
  GenomicloccGreen     220.600     93.202   2.367   0.0179 *  
  GenomiclocdRed       345.383     38.050   9.077  < 2e-16 ***
  GenomicloceOrange    514.300     93.202   5.518 3.47e-08 ***
  GenomiclocfaBlack    360.654     14.115  25.551  < 2e-16 ***
  GenomiclocfBlue      382.567     53.810   7.110 1.20e-12 ***
  GenomiclocgPurple    681.962     57.802  11.798  < 2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#

y3 <- lm(intronnr~Genomicloc, data=intron_total)
y4 <- lm(intronnr~Genomicloc-1, data=intron_total)
anova(y3,y4)
summary(y4)

##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
GenomiclocaAutosome  3.64582    0.02031  179.49  < 2e-16 ***
  GenomiclocbPAR       3.21382    0.16488   19.49  < 2e-16 ***
  GenomicloccGreen     1.80000    0.90910    1.98   0.0477 *  
  GenomiclocdRed       3.71667    0.37114   10.01  < 2e-16 ***
  GenomicloceOrange    4.10000    0.90910    4.51 6.52e-06 ***
  GenomiclocfaBlack    3.91284    0.13768   28.42  < 2e-16 ***
  GenomiclocfBlue      3.90000    0.52487    7.43 1.12e-13 ***
  GenomiclocgPurple    6.65385    0.56380   11.80  < 2e-16 ***

#
y5 <- lm(intronmean~Genomicloc, data=intron_total)
y6 <- lm(intronmean~Genomicloc-1, data=intron_total)
anova(y5,y6)
summary(y6)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
GenomiclocaAutosome   88.360      0.400 220.895  < 2e-16 ***
  GenomiclocbPAR        82.839      3.247  25.512  < 2e-16 ***
  GenomicloccGreen     117.600     17.903   6.569 5.20e-11 ***
  GenomiclocdRed        97.505      7.309  13.340  < 2e-16 ***
  GenomicloceOrange    100.417     17.903   5.609 2.06e-08 ***
  GenomiclocfaBlack     84.171      2.711  31.043  < 2e-16 ***
  GenomiclocfBlue      103.445     10.336  10.008  < 2e-16 ***
  GenomiclocgPurple    107.406     11.103   9.673  < 2e-16 ***
##

#### intron differences
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_intronnrdiffnormalized_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronnrovergenelength)) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.007) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron number difference between A1 and A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_introntotallengdiff_normalized_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronovergenelength)) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.8) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of average intron length between A1 and A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###DE genes and intron size, number
intron_de <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/water_haploid_a1bias_introninfor_geneexpression_v2.txt', header = T)
str(intron_de)

ggplot(intron_de, aes(x=Genomicloc, y=intronnra1, fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,10) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of average intron length between A1 and A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

#logFC and intron difference between a1 and a2.
intron_de_logfc <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/haploidwater_DEnonDE_intron_genomiccompart_fi.txt', header = T)
str(intron_de_logfc)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr, y=logFC.A1.A2, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-6,6) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_intronnrdiff_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr, y=intronnrdiff, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,2) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of intron number between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_intronmeanlength_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr, y=intronmeandiff, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-50,25) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of intron mean length between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_totallength_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr, y=introntotaldiff, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-150,200) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of intron mean length between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###lm stats
y_nr <- lm(logFC.A1.A2 ~ (intronnrdiff + introntotaldiff + intronmeandiff)*chr -1, data = intron_de_logfc) 
summary (y_nr)

#####
Call:
  lm(formula = logFC.A1.A2 ~ (intronnrdiff + introntotaldiff + 
                                intronmeandiff) * chr - 1, data = intron_de_logfc)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.5803  -0.1401   0.0349   0.1896   7.3518 

Coefficients: (2 not defined because of singularities)
Estimate Std. Error t value Pr(>|t|)    
intronnrdiff               -1.442e-01  1.106e-01  -1.303   0.1925    
introntotaldiff            -2.948e-04  9.784e-04  -0.301   0.7632    
intronmeandiff              3.597e-05  1.440e-03   0.025   0.9801    
chraAutosome               -3.433e-02  5.924e-03  -5.795 7.07e-09 ***
  chrbPAR                     5.724e-02  5.008e-02   1.143   0.2530    
chrcGreen                   1.388e-02  3.164e-01   0.044   0.9650    
chrdRed                     2.261e-01  1.176e-01   1.922   0.0546 .  
chreOrange                  6.791e-01  3.211e-01   2.115   0.0344 *  
  chrfaBlack                  5.987e-02  4.409e-02   1.358   0.1745    
chrfBlue                   -6.985e-01  1.541e-01  -4.532 5.91e-06 ***
  chrgPurple                  2.649e-01  1.926e-01   1.375   0.1691    
intronnrdiff:chrbPAR       -2.254e+00  6.782e+00  -0.332   0.7396    
intronnrdiff:chrcGreen             NA         NA      NA       NA    
intronnrdiff:chrdRed        1.486e+00  3.347e+00   0.444   0.6571    
intronnrdiff:chreOrange    -4.653e-01  1.975e+00  -0.236   0.8138    
intronnrdiff:chrfaBlack     3.519e-01  1.713e-01   2.054   0.0400 *  
  intronnrdiff:chrfBlue      -5.220e+00  2.502e+00  -2.086   0.0370 *  
  intronnrdiff:chrgPurple     6.747e-01  5.322e-01   1.268   0.2049    
introntotaldiff:chrbPAR     2.128e-02  7.112e-02   0.299   0.7648    
introntotaldiff:chrcGreen   2.781e-02  1.378e-02   2.018   0.0436 *  
  introntotaldiff:chrdRed    -1.300e-02  3.686e-02  -0.353   0.7243    
introntotaldiff:chreOrange  3.733e-02  1.293e-02   2.887   0.0039 ** 
  introntotaldiff:chrfaBlack -5.468e-05  1.454e-03  -0.038   0.9700    
introntotaldiff:chrfBlue    6.466e-02  3.067e-02   2.108   0.0350 *  
  introntotaldiff:chrgPurple -2.535e-03  5.540e-03  -0.458   0.6473    
intronmeandiff:chrbPAR     -8.958e-02  2.843e-01  -0.315   0.7527    
intronmeandiff:chrcGreen           NA         NA      NA       NA    
intronmeandiff:chrdRed      1.168e-01  1.833e-01   0.637   0.5239    
intronmeandiff:chreOrange  -1.033e-01  4.931e-02  -2.094   0.0362 *  
  intronmeandiff:chrfaBlack   6.006e-03  3.452e-03   1.740   0.0819 .  
intronmeandiff:chrfBlue    -8.706e-02  6.334e-02  -1.374   0.1693    
intronmeandiff:chrgPurple   4.083e-03  1.932e-02   0.211   0.8326    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.5481 on 8874 degrees of freedom
Multiple R-squared:  0.02504,	Adjusted R-squared:  0.02174 
F-statistic: 7.596 on 30 and 8874 DF,  p-value: < 2.2e-16
#####

#expression and intron difference between a1 and a2.
intron_de_expression <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/haploidwater_DEnonDE_intron_genomiccompart_a1a2sep_fi.txt', header = T)
str(intron_de_expression)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_expression_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_expression, aes(x=chr, y=Expression, fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,12) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='LogCPM') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

y_exp <- lm(Expression ~ (intronnr + intronmean + introntotal)*chr -1, data = intron_de_expression) 
summary (y_exp)

#########
lm(formula = Expression ~ (intronnr + intronmean + introntotal) * 
     chr - 1, data = intron_de_expression)
Residuals:
  Min      1Q  Median      3Q     Max 
-9.2190 -1.1548  0.2929  1.3552 10.4056 

Coefficients: (2 not defined because of singularities)
Estimate Std. Error t value Pr(>|t|)    
intronnr                1.598e-01  1.739e-02   9.189  < 2e-16 ***
  intronmean              3.373e-03  4.691e-04   7.189 6.76e-13 ***
  introntotal             6.296e-04  1.879e-04   3.351 0.000806 ***
  chraAutosome            4.133e+00  4.631e-02  89.240  < 2e-16 ***
  chrbPAR                 3.495e+00  3.981e-01   8.780  < 2e-16 ***
  chrcGreen               7.013e+00  1.899e+00   3.694 0.000222 ***
  chrdRed                 3.998e+00  1.098e+00   3.641 0.000272 ***
  chreOrange              2.997e+01  1.199e+01   2.499 0.012480 *  
  chrfaBlack              3.781e+00  3.904e-01   9.684  < 2e-16 ***
  chrfBlue                1.761e+00  3.122e+00   0.564 0.572736    
chrgPurple              5.372e+00  5.245e+00   1.024 0.305729    
intronnr:chrbPAR        2.561e-01  1.787e-01   1.433 0.151921    
intronnr:chrcGreen             NA         NA      NA       NA    
intronnr:chrdRed       -5.093e-02  3.818e-01  -0.133 0.893876    
intronnr:chreOrange    -1.092e+01  4.125e+00  -2.648 0.008100 ** 
  intronnr:chrfaBlack    -1.869e-01  1.438e-01  -1.300 0.193657    
intronnr:chrfBlue      -2.985e-01  1.288e+00  -0.232 0.816683    
intronnr:chrgPurple     4.367e-01  8.431e-01   0.518 0.604436    
intronmean:chrbPAR      3.406e-03  4.555e-03   0.748 0.454613    
intronmean:chrcGreen   -1.109e-02  1.361e-02  -0.815 0.415172    
intronmean:chrdRed      2.308e-03  9.498e-03   0.243 0.808034    
intronmean:chreOrange  -1.980e-01  1.030e-01  -1.923 0.054514 .  
intronmean:chrfaBlack   1.023e-02  4.570e-03   2.239 0.025179 *  
  intronmean:chrfBlue     1.896e-02  3.262e-02   0.581 0.561056    
intronmean:chrgPurple   1.049e-02  4.785e-02   0.219 0.826425    
introntotal:chrbPAR    -2.608e-03  1.940e-03  -1.344 0.178900    
introntotal:chrcGreen          NA         NA      NA       NA    
introntotal:chrdRed     2.878e-03  4.437e-03   0.649 0.516602    
introntotal:chreOrange  8.324e-02  3.155e-02   2.638 0.008351 ** 
  introntotal:chrfaBlack  1.805e-03  1.597e-03   1.130 0.258421    
introntotal:chrfBlue    4.550e-03  1.409e-02   0.323 0.746768    
introntotal:chrgPurple -5.898e-03  8.554e-03  -0.689 0.490540
#########
####

###hemizggous genes and intron######
#####################################
intron_hemi <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/Mvsl_A1_hemi_genomiccompartment.txt', header = T)
str(intron_hemi)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1hemi_intronnr_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi, aes(x=chr, y=intronnr)) + scale_fill_manual(values = c("dark grey")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,10) +              
  theme(legend.position="none") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1hemi_intronmean_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi, aes(x=chr, y=intronmean)) + scale_fill_manual(values = c("dark grey")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,220) +              
  theme(legend.position="none") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron mean length per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1hemi_introntotal_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi, aes(x=chr, y=introntotal)) + scale_fill_manual(values = c("dark grey")) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,800) +              
  theme(legend.position="none") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron total length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###homolog + hemi
intron_hemi_homo <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/Mvsl_A1A2_singleortholog_hemi_genomiccompartment.txt', header = T)
str(intron_hemi_homo)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_homologhemi_intronnr_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi_homo, aes(x=Genomicloc, y=intronnr, fill = ploidy)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("Two copies","Hemizygous"), name="Copy number") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,10) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple","Color")) + 
  labs(x='Genomic compartment', y='Intron total length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_homologhemi_intronmean_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi_homo, aes(x=Genomicloc, y=intronmean, fill = ploidy)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("Two copies","Hemizygous"), name="Copy number") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,210) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple","Color")) + 
  labs(x='Genomic compartment', y='Intron total length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_homologhemi_introntotal_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_hemi_homo, aes(x=Genomicloc, y=introntotal, fill = ploidy)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("Two copies","Hemizygous"), name="Copy number") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1200) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple","Color")) + 
  labs(x='Genomic compartment', y='Intron total length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
