# install packages and load R libraries 
install.packages('ggfortify')
library('ggfortify')

install.packages("missMDA")
library(missMDA)

source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

install.packages("lmerTest")
library(lme4)
library(lmerTest)

#load dataset
PCA_5degen_all <- read.table ("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/PCA/Mvsl_a1a2_exp_compartment_all5traits_nonoriented.txt", header = TRUE)
str(PCA_5degen_all)

PCA_5degen_all1 <- subset(PCA_5degen_all,PCA_5degen_all$youngold!="Centro")
str(PCA_5degen_all1)
#PCA_5degen_all <- data.frame(is.na(PCA_5degen_all))
qqnorm(sqrt(PCA_5degen_all1$abs), pch = 1, frame = FALSE)
qqline(sqrt(PCA_5degen_all1$abs), col = "steelblue")

cor.test(PCA_5degen_all1$abs, PCA_5degen_all1$totalTE, met0hod = "spearm", alternative = "g")

y3 <- glm(sqrt(abs) ~ totalTE + diffGC3 +dn +ratioprot +intronnrdiff+ intronmeandiff + youngold, data =PCA_5degen_all, family =gaussian(link = "identity"))
y4 <- glm(sqrt(abs) ~ totalTE + diffGC3 +dn +ratioprot +intronnrdiff+ intronmeandiff + youngold-1, data =PCA_5degen_all, family =gaussian(link = "identity"))

y3 <- lm(abs ~ diffGC3 +dn +totalTE + ratioprot +intronnrdiff+ intronmeandiff + youngold, data =PCA_5degen_all)
summary (y3)
anova(y3, y4)

y1 <- glm(sqrt(abs) ~ ratioprot*youngold + totalTE*youngold + diffGC3*youngold +dn*youngold +intronnrdiff*youngold+ intronmeandiff*youngold, data = PCA_5degen_all, family=quasi) 
summary(y1)
anova(y1,test="F")
#####
Call:
  Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
TEaverage                            3.067e-03  9.768e-03   0.314 0.753520    
youngoldAuto                         5.375e-01  4.935e-02  10.892  < 2e-16 ***
  youngoldbPAR                         5.977e-02  3.334e-01   0.179 0.857727    
youngoldColorStrata                  1.178e+00  9.985e-01   1.179 0.238287    
youngoldOldStrata                    9.248e-01  3.386e-01   2.731 0.006334 ** 
  averageGC3                          -2.196e-03  8.006e-04  -2.743 0.006102 ** 
  dn                                  -8.178e+00  3.289e+00  -2.486 0.012937 *  
  ds                                   1.221e+01  4.907e+00   2.488 0.012885 *  
  meanprot                            -2.247e-05  7.480e-06  -3.004 0.002676 ** 
  intronmeannr                        -9.023e-03  2.357e-03  -3.828 0.000131 ***
  intronmeantotal                      1.071e-04  2.352e-05   4.553 5.39e-06 ***
  TEaverage:youngoldbPAR              -2.440e-01  1.607e-01  -1.519 0.128911    
TEaverage:youngoldColorStrata        3.321e-01  1.603e-01   2.072 0.038351 *  
  TEaverage:youngoldOldStrata          7.172e-02  4.751e-02   1.510 0.131204    
youngoldbPAR:averageGC3              6.941e-03  5.748e-03   1.207 0.227290    
youngoldColorStrata:averageGC3      -1.678e-02  1.496e-02  -1.122 0.262100    
youngoldOldStrata:averageGC3        -1.497e-03  5.711e-03  -0.262 0.793248    
youngoldbPAR:dn                     -6.191e+01  2.915e+01  -2.123 0.033764 *  
  youngoldColorStrata:dn               1.873e+01  9.881e+00   1.896 0.058016 .  
youngoldOldStrata:dn                 1.037e+01  3.375e+00   3.074 0.002121 ** 
  youngoldbPAR:ds                      1.085e+02  3.539e+01   3.065 0.002185 ** 
  youngoldColorStrata:ds              -3.142e+01  1.614e+01  -1.947 0.051640 .  
youngoldOldStrata:ds                -1.143e+01  5.070e+00  -2.254 0.024222 *  
  youngoldbPAR:meanprot                1.181e-04  7.887e-05   1.497 0.134505    
youngoldColorStrata:meanprot         1.294e-05  9.842e-05   0.131 0.895427    
youngoldOldStrata:meanprot          -3.840e-05  4.928e-05  -0.779 0.435870    
youngoldbPAR:intronmeannr            4.080e-02  2.236e-02   1.825 0.068099 .  
youngoldColorStrata:intronmeannr     1.202e-02  3.692e-02   0.326 0.744719    
youngoldOldStrata:intronmeannr      -2.693e-03  1.481e-02  -0.182 0.855678    
youngoldbPAR:intronmeantotal        -2.919e-04  2.189e-04  -1.333 0.182422    
youngoldColorStrata:intronmeantotal  9.088e-04  2.996e-04   3.033 0.002429 ** 
  youngoldOldStrata:intronmeantotal   -2.450e-04  1.401e-04  -1.749 0.080412 .  
---
#####
y2 <- glm(sqrt(abs) ~ TEaverage + averageGC3 +dn+ ds +meanprot +intronmeannr +intronmeantotal+youngold-1, data = PCA_5degen_all, family=quasi) 
summary(y2)
####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
TEaverage            9.916e-03  9.548e-03   1.039 0.299054    
averageGC3          -2.359e-03  7.879e-04  -2.994 0.002760 ** 
  dn                   1.922e+00  6.913e-01   2.780 0.005454 ** 
  ds                   1.167e+00  1.157e+00   1.009 0.313214    
meanprot            -2.017e-05  7.381e-06  -2.733 0.006289 ** 
  intronmeannr        -8.865e-03  2.321e-03  -3.820 0.000135 ***
  intronmeantotal      1.015e-04  2.310e-05   4.394 1.13e-05 ***
  youngoldAuto         5.463e-01  4.858e-02  11.244  < 2e-16 ***
  youngoldbPAR         5.693e-01  5.238e-02  10.869  < 2e-16 ***
  youngoldColorStrata  6.564e-01  6.355e-02  10.330  < 2e-16 ***
  youngoldOldStrata    7.535e-01  5.304e-02  14.207  < 2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for quasi family taken to be 0.04836812)

Null deviance: 1306.99  on 6114  degrees of freedom
Residual deviance:  295.19  on 6103  degrees of freedom
(2435 observations deleted due to missingness)
AIC: NA
#####


y2 <- glm(sqrt(abs) ~ sqrt(ratioprot):youngold + GC3diff:youngold +dn:youngold+ ds:youngold +intron_nr_diff:youngold +intron_total_diff:youngold +upk5diff+upk10diff:youngold +genediff:youngold, data = PCA_5degen_all, family=gaussian) 
summary(y2)

######
Call:
  glm(formula = sqrt(abs) ~ sqrt(ratioprot):youngold + GC3diff:youngold + 
        dn:youngold + ds:youngold + intron_nr_diff:youngold + intron_total_diff:youngold + 
        upk5diff + upk10diff:youngold + genediff:youngold, family = quasi, 
      data = PCA_5degen_all)

Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-0.84708  -0.14648  -0.02776   0.11469   2.00551  

Coefficients: (3 not defined because of singularities)
Estimate Std. Error t value Pr(>|t|)    
(Intercept)                           -1.117e+00  3.352e-01  -3.331 0.000870 ***
  upk5diff                              -2.496e-03  7.280e-03  -0.343 0.731737    
sqrt(ratioprot):youngoldAuto           1.511e+00  3.352e-01   4.507 6.69e-06 ***
  sqrt(ratioprot):youngoldbPAR           1.535e+00  3.362e-01   4.566 5.07e-06 ***
  sqrt(ratioprot):youngoldColorStrata    1.500e+00  3.395e-01   4.418 1.02e-05 ***
  sqrt(ratioprot):youngoldOldStrata      1.676e+00  3.342e-01   5.014 5.49e-07 ***
  youngoldAuto:GC3diff                  -5.014e-02  3.717e-02  -1.349 0.177375    
youngoldbPAR:GC3diff                   6.974e+00  8.969e+00   0.778 0.436810    
youngoldColorStrata:GC3diff            1.148e-01  7.028e-02   1.633 0.102461    
youngoldOldStrata:GC3diff              6.760e-02  1.402e-02   4.822 1.46e-06 ***
  youngoldAuto:dn                       -7.129e+00  3.519e+00  -2.026 0.042794 *  
  youngoldbPAR:dn                       -2.477e+02  2.257e+02  -1.097 0.272532    
youngoldColorStrata:dn                 3.671e+01  1.026e+01   3.578 0.000349 ***
  youngoldOldStrata:dn                   2.468e+00  7.351e-01   3.357 0.000793 ***
  youngoldAuto:ds                        1.140e+01  5.046e+00   2.259 0.023931 *  
  youngoldbPAR:ds                        1.231e+02  3.603e+01   3.418 0.000634 ***
  youngoldColorStrata:ds                -4.455e+01  1.728e+01  -2.578 0.009970 ** 
  youngoldOldStrata:ds                   2.835e+00  1.337e+00   2.121 0.033973 *  
  youngoldAuto:intron_nr_diff            1.101e-02  5.073e-02   0.217 0.828132    
youngoldbPAR:intron_nr_diff                   NA         NA      NA       NA    
youngoldColorStrata:intron_nr_diff     3.740e-02  2.077e-01   0.180 0.857090    
youngoldOldStrata:intron_nr_diff       7.675e-02  3.317e-02   2.314 0.020711 *  
  youngoldAuto:intron_total_diff        -4.295e-04  3.833e-04  -1.121 0.262522    
youngoldbPAR:intron_total_diff         7.568e-02  7.822e-02   0.967 0.333335    
youngoldColorStrata:intron_total_diff  2.278e-04  8.575e-04   0.266 0.790517    
youngoldOldStrata:intron_total_diff    7.216e-04  2.288e-04   3.154 0.001619 ** 
  youngoldAuto:upk10diff                -8.209e-02  2.628e-02  -3.123 0.001796 ** 
  youngoldbPAR:upk10diff                 2.313e-01  2.198e-01   1.053 0.292525    
youngoldColorStrata:upk10diff         -3.354e-02  3.091e-02  -1.085 0.277926    
youngoldOldStrata:upk10diff           -2.915e-02  9.753e-03  -2.989 0.002812 ** 
  youngoldAuto:genediff                  3.270e-04  1.559e-02   0.021 0.983272    
youngoldbPAR:genediff                         NA         NA      NA       NA    
youngoldColorStrata:genediff                  NA         NA      NA       NA    
youngoldOldStrata:genediff             2.813e-02  4.095e-02   0.687 0.492130    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for quasi family taken to be 0.04764018)

Null deviance: 314.06  on 6113  degrees of freedom
Residual deviance: 289.80  on 6083  degrees of freedom
(2432 observations deleted due to missingness)
AIC: NA

Number of Fisher Scoring iterations: 2
######

##using variables of degeneration traits as column names
PCA_5degen1 <- as.data.frame(scale(PCA_5degen_all[,7:25]))
str(PCA_5degen1)
sum(is.na(PCA_5degen1))

pc <- pca(PCA_5degen1, nPcs=5, method="ppca")
imputed <- completeObs(pc)

####
pca <- prcomp(imputed, center = TRUE, scale. =  TRUE)
summary(pca)

pca$gencomp <- PCA_5degen_all[,2]
pca$de2 <- PCA_5degen_all[,6]
str(pca)

ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,  labels=rownames(mtcars), groups=mtcars.country) +
  scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue"))+
  ggtitle("PCA of mtcars dataset")+
  theme_minimal()+
  theme(legend.position = "bottom")

#pca1 and pca2
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca1pca2.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

#pca3 and pca4
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca3pca4.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, choices=c(3,4), obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

#pca5 and pca6
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca5pca6.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, choices=c(5,6), obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

