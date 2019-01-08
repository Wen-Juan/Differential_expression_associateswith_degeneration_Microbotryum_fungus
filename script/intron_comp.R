#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
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

###scatter plot on positive correlation between intron length and coding sequencing length.
ggplot(intron_total, aes(x=introntotal, y=genelength)) + geom_point()

#plot intron mean length
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_intronnr_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronnr, fill=Haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Average intron number per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

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
intron_de <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/haploidwater_DEnonDE_intron_genomiccompart_a1a2sep_fi.txt', header = T)
str(intron_de)

ggplot(intron_de, aes(x=chr, y=introntotal, fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1500) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Difference of average intron length between A1 and A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

#logFC and intron difference between a1 and a2.
intron_de_logfc <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/waterhaploid_Logcpm_0.05_DEnonDE_exp_fi.txt', header = T)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr1, y=logFC.A1.A2, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey","dodgerblue2",), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-6,6) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_intronnrdiffnormalize_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr1, y=introntotaldiff, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-150,150) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference of normalized intron number between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_intronmeanlength_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr1, y=intronmeandiff, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.3,0.6) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference of intron mean length between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_DEnonDE_totallengthnormalize_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_de_logfc, aes(x=chr1, y=introntotallengthnormalize, fill=DE_status)) + scale_fill_manual(values = c("firebrick2","grey"), labels=c("DE","Non-DE"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.75) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference of intron mean length between a1 and a2') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
####
NRR <- subset(intron_de_logfc, intron_de_logfc$chr1 == "NRR")
wilcox.test(NRR$introntotallengthnormalize[NRR$DE_status=='DE'],NRR$introntotallengthnormalize[NRR$DE_status=='NON'],exact = FALSE)
#W = 4944, p-value = 0.3087
wilcox.test(NRR$intronnrnormalize[NRR$DE_status=='DE'],NRR$intronnrnormalize[NRR$DE_status=='NON'],exact = FALSE)
#W = 4688.5, p-value = 0.6991

###lm stats
y_nr <- lm(logFC.A1.A2 ~ (intronnrnormalize + introntotallengthnormalize + intronmeandiff)*chr -1, data = intron_de_logfc) 
summary (y_nr)

#####
lm(formula = logFC.A1.A2 ~ (intronnrnormalize + introntotallengthnormalize + 
                              intronmeandiff) * chr - 1, data = intron_de_logfc)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.3875  -0.1449   0.0339   0.1911   7.3548 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
intronnrnormalize                     -8.038e+00  5.602e+00  -1.435 0.151329    
introntotallengthnormalize             2.279e-01  5.149e-02   4.426 9.70e-06 ***
  intronmeandiff                         2.674e-04  1.184e-03   0.226 0.821254    
chraAutosome                          -6.238e-02  1.237e-02  -5.044 4.65e-07 ***
  chrbPAR                               -4.079e-03  1.090e-01  -0.037 0.970154    
chrcGreen                              7.729e-02  9.071e-01   0.085 0.932095    
chrdRed                                5.472e-02  2.607e-01   0.210 0.833748    
chreOrange                            -5.097e+00  1.072e+00  -4.753 2.03e-06 ***
  chrfaBlack                            -7.135e-03  8.557e-02  -0.083 0.933551    
chrfBlue                               1.323e+00  4.039e-01   3.276 0.001059 ** 
  chrgPurple                             1.912e-01  5.137e-01   0.372 0.709715    
intronnrnormalize:chrbPAR             -1.467e+01  4.082e+01  -0.359 0.719237    
intronnrnormalize:chrcGreen            3.500e+02  5.785e+02   0.605 0.545181    
intronnrnormalize:chrdRed              1.210e+02  1.001e+02   1.209 0.226680    
intronnrnormalize:chreOrange          -1.038e+04  2.713e+03  -3.827 0.000130 ***
  intronnrnormalize:chrfaBlack          -6.633e+01  4.198e+01  -1.580 0.114095    
intronnrnormalize:chrfBlue            -1.325e+03  2.166e+02  -6.117 9.91e-10 ***
  intronnrnormalize:chrgPurple          -5.175e+02  4.316e+02  -1.199 0.230549    
introntotallengthnormalize:chrbPAR     3.753e-01  3.498e-01   1.073 0.283315    
introntotallengthnormalize:chrcGreen  -4.283e+00  9.981e+00  -0.429 0.667832    
introntotallengthnormalize:chrdRed    -7.996e-01  7.345e-01  -1.089 0.276385    
introntotallengthnormalize:chreOrange  1.083e+02  2.090e+01   5.179 2.28e-07 ***
  introntotallengthnormalize:chrfaBlack  1.010e+00  4.046e-01   2.497 0.012537 *  
  introntotallengthnormalize:chrfBlue    4.797e+00  2.146e+00   2.235 0.025430 *  
  introntotallengthnormalize:chrgPurple  4.623e+00  3.719e+00   1.243 0.213919    
intronmeandiff:chrbPAR                -5.741e-03  3.989e-03  -1.439 0.150192    
intronmeandiff:chrcGreen               2.135e-02  1.277e-01   0.167 0.867237    
intronmeandiff:chrdRed                 4.099e-02  1.032e-01   0.397 0.691255    
intronmeandiff:chreOrange              2.645e-01  3.566e-02   7.418 1.30e-13 ***
  intronmeandiff:chrfaBlack              5.467e-03  2.294e-03   2.384 0.017165 *  
  intronmeandiff:chrfBlue                3.014e-02  8.811e-03   3.421 0.000627 ***
  intronmeandiff:chrgPurple             -1.534e-02  8.754e-03  -1.752 0.079820 .  
#####
y_nr <- lm(logFC.A1.A2 ~ (intronnrnormalize + introntotallengthnormalize + intronmeandiff)*chr1 -1, data = intron_de_logfc) 
summary (y_nr)
lm(formula = logFC.A1.A2 ~ (intronnrnormalize + introntotallengthnormalize + 
                              intronmeandiff) * chr1 - 1, data = intron_de_logfc)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.4618  -0.1451   0.0339   0.1905   7.3548 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
intronnrnormalize                   -8.038e+00  5.660e+00  -1.420   0.1556    
introntotallengthnormalize       2.279e-01  5.203e-02   4.381 1.20e-05 ***
  intronmeandiff                       2.674e-04  1.196e-03   0.224   0.8231    
chr1aAutosome                       -6.238e-02  1.250e-02  -4.992 6.08e-07 ***
  chr1bPAR                            -4.079e-03  1.101e-01  -0.037   0.9705    
chr1NRR                              9.323e-02  7.735e-02   1.205   0.2281    
intronnrnormalize:chr1bPAR          -1.467e+01  4.124e+01  -0.356   0.7220    
intronnrnormalize:chr1NRR           -8.174e+01  3.719e+01  -2.198   0.0280 *  
  introntotallengthnormalize:chr1bPAR  3.753e-01  3.534e-01   1.062   0.2883    
introntotallengthnormalize:chr1NRR   7.222e-01  3.402e-01   2.123   0.0338 *  
  intronmeandiff:chr1bPAR             -5.741e-03  4.031e-03  -1.424   0.1544    
intronmeandiff:chr1NRR               6.616e-03  2.213e-03   2.989   0.0028 ** 
  
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

y1_nr <- lm(logFC.A1.A2 ~ intronmeandiff*chr -1, data = intron_de_logfc) 
summary (y1_nr)

####
lm(formula = logFC.A1.A2 ~ intronmeandiff * chr - 1, data = intron_de_logfc)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.5808  -0.1404   0.0347   0.1896   7.3518 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
intronmeandiff             0.000300   0.001192   0.252   0.8013    
chraAutosome              -0.034298   0.005944  -5.770 8.20e-09 ***
  chrbPAR                    0.062265   0.049785   1.251   0.2111    
chrcGreen                  0.013879   0.317609   0.044   0.9651    
chrdRed                    0.209264   0.107897   1.939   0.0525 .  
chreOrange                 1.769476   0.281762   6.280 3.54e-10 ***
  chrfaBlack                 0.060356   0.043970   1.373   0.1699    
chrfBlue                  -0.622661   0.150475  -4.138 3.54e-05 ***
  chrgPurple                 0.211819   0.174990   1.210   0.2261    
intronmeandiff:chrbPAR    -0.005661   0.004000  -1.415   0.1571    
intronmeandiff:chrcGreen   0.054776   0.027644   1.981   0.0476 *  
  intronmeandiff:chrdRed     0.046258   0.103099   0.449   0.6537    
intronmeandiff:chreOrange  0.036909   0.019853   1.859   0.0630 .  
intronmeandiff:chrfaBlack  0.005871   0.002294   2.559   0.0105 *  
  intronmeandiff:chrfBlue    0.044433   0.008620   5.155 2.60e-07 ***
  intronmeandiff:chrgPurple -0.011486   0.008215  -1.398   0.1621 
#######
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
ggplot(intron_hemi_homo, aes(x=Genomicloc, y=intronnr, fill = Haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("Two copies","Hemizygous"), name="Copy number") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,10) +              
  #theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple","Color")) + 
  labs(x='Genomic compartment', y='Average intron number per gene') +
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
