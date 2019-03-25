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
PCA_5degen_all <- read.table ("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/PCA/Mvsl_all5traits_withmissingvalues_oriented_all_fi.txt", header = TRUE)
str(PCA_5degen_all)
qqnorm(sqrt(PCA_5degen_all$abs), pch = 1, frame = FALSE)
qqline(sqrt(PCA_5degen_all$abs), col = "steelblue")

y1 <- glm(sqrt(abs) ~ sqrt(ratioprot)*youngold + GC3diff*youngold +dn*youngold+ ds*youngold +intron_nr_diff*youngold +intron_total_diff*youngold +genetoup10kdiff*youngold-1, data = PCA_5degen_all, family=quasi) 
summary(y1)
#####
Call:
  glm(formula = sqrt(abs) ~ sqrt(ratioprot) * youngold + GC3diff * 
        youngold + dn * youngold + ds * youngold + intron_nr_diff * 
        youngold + intron_total_diff * youngold + genetoup10kdiff * 
        youngold - 1, family = quasi, data = PCA_5degen_all)

Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-0.79624  -0.14657  -0.02725   0.11452   2.00541  

Coefficients: (1 not defined because of singularities)
Estimate Std. Error t value Pr(>|t|)    
sqrt(ratioprot)                       -4.022e-01  6.459e-01  -0.623  0.53348    
youngoldAuto                           7.964e-01  6.459e-01   1.233  0.21764    
youngoldbPAR                          -2.676e+10  4.173e+10  -0.641  0.52132    
youngoldColorStrata                   -4.713e+00  4.087e+00  -1.153  0.24889    
youngoldOldStrata                     -2.288e+00  3.550e-01  -6.444 1.25e-10 ***
  GC3diff                               -6.841e-02  3.021e-02  -2.264  0.02359 *  
  dn                                    -7.450e+00  3.309e+00  -2.251  0.02440 *  
  ds                                     1.169e+01  4.967e+00   2.354  0.01858 *  
  intron_nr_diff                         4.905e-02  3.913e-02   1.253  0.21015    
intron_total_diff                     -4.620e-05  2.488e-04  -0.186  0.85270    
genetoup10kdiff                       -1.782e-02  1.200e-02  -1.485  0.13761    
sqrt(ratioprot):youngoldbPAR           2.676e+10  4.173e+10   0.641  0.52132    
sqrt(ratioprot):youngoldColorStrata    5.547e+00  4.121e+00   1.346  0.17832    
sqrt(ratioprot):youngoldOldStrata      3.252e+00  7.355e-01   4.422 9.94e-06 ***
  youngoldbPAR:GC3diff                   1.911e+09  2.980e+09   0.641  0.52132    
youngoldColorStrata:GC3diff            3.151e-01  1.125e-01   2.801  0.00511 ** 
  youngoldOldStrata:GC3diff              1.003e-01  3.334e-02   3.008  0.00264 ** 
  youngoldbPAR:dn                       -4.778e+10  7.450e+10  -0.641  0.52132    
youngoldColorStrata:dn                 2.042e+01  9.721e+00   2.101  0.03572 *  
  youngoldOldStrata:dn                   9.601e+00  3.387e+00   2.835  0.00460 ** 
  youngoldbPAR:ds                        1.114e+02  3.622e+01   3.076  0.00211 ** 
  youngoldColorStrata:ds                -4.702e+00  1.632e+01  -0.288  0.77330    
youngoldOldStrata:ds                  -7.640e+00  5.133e+00  -1.489  0.13664    
youngoldbPAR:intron_nr_diff                   NA         NA      NA       NA    
youngoldColorStrata:intron_nr_diff     5.902e-01  1.260e-01   4.686 2.85e-06 ***
  youngoldOldStrata:intron_nr_diff       1.576e-03  4.704e-02   0.033  0.97328    
youngoldbPAR:intron_total_diff         1.200e-01  8.068e-02   1.487  0.13694    
youngoldColorStrata:intron_total_diff  3.053e-03  6.319e-04   4.831 1.39e-06 ***
  youngoldOldStrata:intron_total_diff    8.356e-04  3.049e-04   2.740  0.00616 ** 
  youngoldbPAR:genetoup10kdiff           2.359e-01  1.563e-01   1.509  0.13125    
youngoldColorStrata:genetoup10kdiff    3.159e-02  2.744e-02   1.151  0.24972    
youngoldOldStrata:genetoup10kdiff      4.587e-03  1.304e-02   0.352  0.72492    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for quasi family taken to be 0.04723679)

Null deviance: 1306.99  on 6114  degrees of freedom
Residual deviance:  287.34  on 6083  degrees of freedom
(2432 observations deleted due to missingness)
AIC: NA

Number of Fisher Scoring iterations: 2
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

