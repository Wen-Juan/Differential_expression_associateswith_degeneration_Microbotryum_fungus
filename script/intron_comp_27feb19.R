#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load data to reandomdize the non-DE gene expression direction, on 15thFeb.2019
nonDE_intron <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/intron_random_diff.txt', header = T)
str(nonDE_intron)

intron_nr <- data.frame(nonDE_intron$intron_nr_a1, nonDE_intron$intron_nr_a2)
intron_mean <- data.frame(nonDE_intron$intron_mean_a1, nonDE_intron$intron_mean_a2)
intron_total <- data.frame(nonDE_intron$intron_total_a1, nonDE_intron$intron_total_a2)

intron_nr_rand <- randomizeMatrix(intron_nr,null.model = "richness",iterations = 1000)
intron_mean_rand <- randomizeMatrix(intron_mean,null.model = "richness",iterations = 1000)
intron_total_rand <- randomizeMatrix(intron_total,null.model = "richness",iterations = 1000)

total1 <- cbind(intron_nr_rand,intron_mean_rand, intron_total_rand)
head(total1)

total_data <- cbind(nonDE_intron, total1)
write.table(total_data,file = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_compart_randomdized.txt",quote=F, row.names=T, sep='\t')

#intron_random <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_exp_compart_randomdize.txt', header = T)
#str(intron_random)

intron_random <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_exp_compart_low-high.txt', header = T)
str(intron_random)

intron_random_rmcentro <- subset(intron_random, intron_random$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_introntotal_diff1_youngold.pdf", width=8, height=8)
pa <- ggplot(intron_random_rmcentro, aes(x=youngold, y=intron_total_diff, fill=DE2)) + 
  scale_fill_manual(values = c("white","grey"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-250,250) +   
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length difference') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

intron_random_rmcentro_old <- subset(intron_random_rmcentro, intron_random_rmcentro$youngold == "OldStrata")
wilcox.test(intron_random_rmcentro_old$intron_mean_diff[intron_random_rmcentro_old$DE2=="DE"], intron_random_rmcentro_old$intron_mean_diff[intron_random_rmcentro_old$DE2=="NON"],exact=FALSE)
#W = 4243, p-value = 0.9264 
wilcox.test(intron_random_rmcentro_old$intron_total_diff[intron_random_rmcentro_old$DE2=="DE"], intron_random_rmcentro_old$intron_total_diff[intron_random_rmcentro_old$DE2=="NON"], exact=FALSE)
#W = 4694.5, p-value = 0.1813
wilcox.test(intron_random_rmcentro_old$intron_nr_diff[intron_random_rmcentro_old$DE2=="DE"], intron_random_rmcentro_old$intron_nr_diff[intron_random_rmcentro_old$DE2=="NON"], exact=FALSE)
#W = 4325, p-value = 0.6985

intron_random_rmcentro_young <- subset(intron_random_rmcentro, intron_random_rmcentro$youngold == "ColorStrata")
wilcox.test(intron_random_rmcentro_young$intron_mean_diff[intron_random_rmcentro_young$DE2=="DE"], intron_random_rmcentro_young$intron_mean_diff[intron_random_rmcentro_young$DE2=="NON"], exact=FALSE)
#W = 82, p-value = 0.007993
wilcox.test(intron_random_rmcentro_young$intron_total_diff[intron_random_rmcentro_young$DE2=="DE"], intron_random_rmcentro_young$intron_total_diff[intron_random_rmcentro_young$DE2=="NON"], exact=FALSE)
#W = 84.5, p-value = 0.004592
wilcox.test(intron_random_rmcentro_young$intron_nr_diff[intron_random_rmcentro_young$DE2=="DE"], intron_random_rmcentro_young$intron_nr_diff[intron_random_rmcentro_young$DE2=="NON"], exact=FALSE)
#W = 60.5, p-value = 0.09752

intron_random_rmcentro_auto <- subset(intron_random_rmcentro, intron_random_rmcentro$youngold == "Auto")
intron_random_rmcentro_par <- subset(intron_random_rmcentro, intron_random_rmcentro$youngold == "bPAR")
wilcox.test(intron_random_rmcentro_auto$intron_mean_diff[intron_random_rmcentro_auto$DE2=="DE"], intron_random_rmcentro_auto$intron_mean_diff[intron_random_rmcentro_auto$DE2=="NON"], exact=FALSE)
#W = 1946300, p-value = 0.7075
wilcox.test(intron_random_rmcentro_auto$intron_total_diff[intron_random_rmcentro_auto$DE2=="DE"], intron_random_rmcentro_auto$intron_total_diff[intron_random_rmcentro_auto$DE2=="NON"], exact=FALSE)
#W = 1983900, p-value = 0.03227
wilcox.test(intron_random_rmcentro_auto$intron_nr_diff[intron_random_rmcentro_auto$DE2=="DE"], intron_random_rmcentro_auto$intron_nr_diff[intron_random_rmcentro_auto$DE2=="NON"], exact=FALSE)
#W = 1982000, p-value = 1.979e-05

wilcox.test(intron_random_rmcentro_par$intron_mean_diff[intron_random_rmcentro_par$DE2=="DE"], intron_random_rmcentro_par$intron_mean_diff[intron_random_rmcentro_par$DE2=="NON"], exact=FALSE)
#W = 677, p-value = 0.2248
wilcox.test(intron_random_rmcentro_par$intron_total_diff[intron_random_rmcentro_par$DE2=="DE"], intron_random_rmcentro_par$intron_total_diff[intron_random_rmcentro_par$DE2=="NON"], exact=FALSE)
#W = 677.5, p-value = 0.2212
wilcox.test(intron_random_rmcentro_par$intron_nr_diff[intron_random_rmcentro_par$DE2=="DE"], intron_random_rmcentro_par$intron_nr_diff[intron_random_rmcentro_par$DE2=="NON"], exact=FALSE)
#W = 612, p-value = NA


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronmean_diff1_youngold.pdf", width=8, height=8)
pb <- ggplot(intron_random_rmcentro, aes(x=youngold, y=intron_mean_diff, fill=DE2)) + 
  scale_fill_manual(values = c("white","grey"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-20,20) +       
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Mean intron length difference') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronnumber_diff1_youngold.pdf", width=8, height=8)
pc <- ggplot(intron_random_rmcentro, aes(x=youngold, y=intron_nr_diff, fill=DE2)) + 
  scale_fill_manual(values = c("white","grey"), labels=c("DE","Non-DE"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +             
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Intron number difference') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intron_combine_3figs_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(pa,pb,pc,labels=c('A','B','C'))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_intronnr_corr_youngold.pdf", width=8, height=8)
pa1 <- ggplot(intron_random_rmcentro, aes(x=intron_nr_diff, y=abs,color=DE2,shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"),guide = FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  xlim(-8,8) +
  ylim(0,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='Intron number difference', y= 'Gene expression ratio |Log2(A1/A2)|')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_intronmeanlength_corr_youngold.pdf", width=8, height=8)
pa2 <- ggplot(intron_random_rmcentro, aes(x=intron_mean_diff, y=abs,color=DE2,shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  xlim(-200,200) +
  ylim(0,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='Intron mean length difference', y= 'Gene expression ratio |Log2(A1/A2)|')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_introntotallength_corr_youngold.pdf", width=8, height=8)
pa3 <-ggplot(intron_random_rmcentro, aes(x=intron_total_diff, y=abs,color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  xlim(-600,600) +
  ylim(0,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='Intron total length difference', y= 'Gene expression ratio |Log2(A1/A2)|')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intron_3figures.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(pa1, pa3, pa2,labels=c('A','B','C'))
dev.off()

cor.test(intron_random_rmcentro$intron_nr_diff[intron_random_rmcentro$DE2 == "NON"], intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "NON"], method=c("pearson"))
cor.test(intron_random_rmcentro$intron_nr_diff[intron_random_rmcentro$DE2 == "DE"], intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"], method=c("pearson"))
#Pearson's product-moment correlation

data:  intron_random_rmcentro$intron_nr_diff[intron_random_rmcentro$DE2 ==  and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]    "DE"] and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]
t = -1.3069, df = 589, p-value = 0.1918
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  -0.13384218  0.02699846
sample estimates:
  cor 
-0.05377062 
###

cor.test(intron_random_rmcentro$intron_mean_diff[intron_random_rmcentro$DE2 == "DE"], intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"], method=c("pearson"))
##Pearson's product-moment correlation

data:  intron_random_rmcentro$intron_mean_diff[intron_random_rmcentro$DE2 ==  and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]    "DE"] and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]
t = -0.17708, df = 589, p-value = 0.8595
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  -0.08789672  0.07339873
sample estimates:
  cor 
-0.007296456 
###

cor.test(intron_random_rmcentro$intron_total_diff[intron_random_rmcentro$DE2 == "DE"], intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"], method=c("pearson"))
###Pearson's product-moment correlation
data:  intron_random_rmcentro$intron_total_diff[intron_random_rmcentro$DE2 ==  and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]    "DE"] and intron_random_rmcentro$abs[intron_random_rmcentro$DE2 == "DE"]
t = -0.78764, df = 589, p-value = 0.4312
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  -0.1127938  0.0483416
sample estimates:
  cor 
-0.03243686
###


y1a <- lm(abs ~ DE2*intron_nr_diff-1, data=intron_random_rmcentro)
summary(y1a)
###
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE                  1.889194   0.018142 104.132  < 2e-16 ***
  DE2NON                 0.214065   0.004943  43.303  < 2e-16 ***
  intron_nr_diff        -0.191351   0.043042  -4.446 8.87e-06 ***
  DE2NON:intron_nr_diff  0.201402   0.053830   3.741 0.000184 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4409 on 8541 degrees of freedom
Multiple R-squared:  0.5983,	Adjusted R-squared:  0.5981 
F-statistic:  3180 on 4 and 8541 DF,  p-value: < 2.2e-16
###

y1 <- lm(abs ~ DE2/intron_nr_diff-1, data=intron_random_rmcentro)
summary(y1)
###
Residuals:
  Min      1Q  Median      3Q     Max 
-1.3549 -0.1541 -0.0676  0.0819  9.8226 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE                  1.889194   0.018142 104.132  < 2e-16 ***
  DE2NON                 0.214065   0.004943  43.303  < 2e-16 ***
  DE2DE:intron_nr_diff  -0.191351   0.043042  -4.446 8.87e-06 ***
  DE2NON:intron_nr_diff  0.010050   0.032327   0.311    0.756    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4409 on 8541 degrees of freedom
Multiple R-squared:  0.5983,	Adjusted R-squared:  0.5981 
F-statistic:  3180 on 4 and 8541 DF,  p-value: < 2.2e-16 
###

y2a <- lm(abs ~ DE2*intron_mean_diff-1, data=intron_random_rmcentro)
summary(y2a)
###
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE                    1.8870913  0.0181578 103.927   <2e-16 ***
  DE2NON                   0.2140502  0.0049490  43.251   <2e-16 ***
  intron_mean_diff        -0.0009741  0.0016165  -0.603    0.547    
DE2NON:intron_mean_diff  0.0012432  0.0018356   0.677    0.498    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4414 on 8541 degrees of freedom
Multiple R-squared:  0.5974,	Adjusted R-squared:  0.5972 
F-statistic:  3168 on 4 and 8541 DF,  p-value: < 2.2e-16
###

y2 <- lm(abs ~ DE2/intron_mean_diff-1, data=intron_random_rmcentro)
summary(y2)
###
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE                    1.8870913  0.0181579 103.927   <2e-16 ***
  DE2NON                   0.2140603  0.0049490  43.253   <2e-16 ***
  DE2DE:intron_mean_diff  -0.0009741  0.0016165  -0.603    0.547    
DE2NON:intron_mean_diff -0.0001214  0.0008698  -0.140    0.889    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4414 on 8541 degrees of freedom
Multiple R-squared:  0.5974,	Adjusted R-squared:  0.5972 
F-statistic:  3168 on 4 and 8541 DF,  p-value: < 2.2e-16
####

y3 <- lm(abs ~ DE2*intron_total_diff-1, data=intron_random_rmcentro)
summary(y3)
####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE                     1.8893463  0.0181710 103.976  < 2e-16 ***
  DE2NON                    0.2140534  0.0049470  43.269  < 2e-16 ***
  intron_total_diff        -0.0010153  0.0003789  -2.680  0.00738 ** 
  DE2NON:intron_total_diff  0.0009714  0.0004221   2.301  0.02140 *  
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4412 on 8541 degrees of freedom
Multiple R-squared:  0.5977,	Adjusted R-squared:  0.5975 
F-statistic:  3172 on 4 and 8541 DF,  p-value: < 2.2e-16
####



####15feb.2019, check coding sequence length and protein lengh variation across genomic compartment
cds_prot <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_cds_prot_compart.txt', header = T)
str(cds_prot)

cds_prot_rmcentro <- subset(cds_prot, cds_prot$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_cds_combined_youngold.pdf", width=8, height=8)
pa1 <- ggplot(cds_prot_rmcentro, aes(x=youngold, y=cds, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"), guide = FALSE) + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,4000) +                 
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Coding sequence length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


cds_prot_rmcentro_auto <- subset(cds_prot_rmcentro, cds_prot_rmcentro$youngold =="Auto")
cds_prot_rmcentro_old <- subset(cds_prot_rmcentro, cds_prot_rmcentro$youngold =="OldStrata")
wilcox.test(cds_prot_rmcentro_auto$cds[cds_prot_rmcentro_auto$DE2=="DE"], cds_prot_rmcentro_auto$cds[cds_prot_rmcentro_auto$DE2=="NON"],extact=FALSE)
#W = 7541200, p-value = 0.06819
wilcox.test(cds_prot_rmcentro_auto$prot[cds_prot_rmcentro_auto$DE2=="DE"], cds_prot_rmcentro_auto$prot[cds_prot_rmcentro_auto$DE2=="NON"],extact=FALSE)
#W = 7541300, p-value = 0.06819
wilcox.test(cds_prot_rmcentro_old$cds[cds_prot_rmcentro_old$DE2=="DE"], cds_prot_rmcentro_old$cds[cds_prot_rmcentro_old$DE2=="NON"],extact=FALSE)
#W = 16276, p-value = 0.5894
wilcox.test(cds_prot_rmcentro_old$prot[cds_prot_rmcentro_old$DE2=="DE"], cds_prot_rmcentro_old$prot[cds_prot_rmcentro_old$DE2=="NON"],extact=FALSE)
#W = 16276, p-value = 0.5894

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_prot_combined_youngold.pdf", width=8, height=8)
pa2 <- ggplot(cds_prot_rmcentro, aes(x=youngold, y=prot, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"), labels=c("DE","Non-DE"), name="Bias direction") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1800) +             
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Protein length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_coding_proteinratio_combine_figs.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(pa1, pa2, labels=c('A','B'))
dev.off()


cds_prot_divide <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_cds_prot_divide_compart2.txt', header = T)
str(cds_prot_divide)

cds_prot_divide_rmcentro <- subset(cds_prot_divide, cds_prot_divide$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronnr_normalize_youngold.pdf", width=8, height=8)
ggplot(cds_prot_divide_rmcentro, aes(x=youngold, y=intron_nr_byprota1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","grey"), labels=c("DE","Non-DE"), name="Bias Direction") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.03) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Intron number normalized by protein length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

cds_prot_divide_auto <- subset(cds_prot_divide, cds_prot_divide$youngold=="Auto")
cds_prot_divide_old <- subset(cds_prot_divide, cds_prot_divide$youngold=="OldStrata")
wilcox.test(cds_prot_divide_auto$intron_nr_byprota1[cds_prot_divide_auto$DE2=="DE"], cds_prot_divide_auto$intron_nr_byprota1[cds_prot_divide_auto$DE2=="NON"], exact = FALSE)
#W = 7968900, p-value = 0.2704
wilcox.test(cds_prot_divide_auto$intron_total_byprota1[cds_prot_divide_auto$DE2=="DE"], cds_prot_divide_auto$intron_total_byprota1[cds_prot_divide_auto$DE2=="NON"], exact = FALSE)
#W = 7961900, p-value = 0.2916

wilcox.test(cds_prot_divide_old$intron_nr_byprota1[cds_prot_divide_old$DE2=="DE"], cds_prot_divide_old$intron_nr_byprota1[cds_prot_divide_old$DE2=="NON"], exact = FALSE)
#W = 15562, p-value = 0.2193
wilcox.test(cds_prot_divide_old$intron_total_byprota1[cds_prot_divide_old$DE2=="DE"], cds_prot_divide_old$intron_total_byprota1[cds_prot_divide_old$DE2=="NON"], exact = FALSE)
#W = 15154, p-value = 0.1047

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_introntotallength_normalize_youngold.pdf", width=8, height=8)
ggplot(cds_prot_divide_rmcentro, aes(x=youngold, y=intron_total_byprota1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","grey"), labels=c("DE","Non-DE"), name="Bias Direction") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,3) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Intron totral length normalized by protein length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


#load the corresponding data files, on 05feb2019
intron_total <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_homolog_intron.txt', header = T)
str(intron_total)
intron_total_rmcentro <- subset(intron_total, intron_total$youngold!="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronnr_diff_youngold.pdf", width=8, height=8)
p1 <- ggplot(intron_total_rmcentro, aes(x=youngold, y=intron_nr_a1, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                 
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronmean_diff_youngold.pdf", width=8, height=8)
p2 <- ggplot(intron_total_rmcentro, aes(x=youngold, y=intron_mean_a1, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"),guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,220) +   
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_introntotallength_diff_youngold.pdf", width=8, height=8)
p3 <- ggplot(intron_total_rmcentro, aes(x=youngold, y=intron_total_a1, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1300) +     
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intron_all_3figs_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(p1,p2,p3,labels=c('A','B','C'))
dev.off()

y <- lm(intron_total_a1~youngold*haploid, data=intron_total_rmcentro)
anova(y)
summary(y)

y1 <- lm(intron_mean_a1~youngold*haploid, data=intron_total_rmcentro)
anova(y1)
summary(y1)

y3 <- lm(intron_nr_a1~youngold*haploid, data=intron_total_rmcentro)
anova(y3)
summary(y3)

intron_exp <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_exp_compart.txt', header = T)
str(intron_exp)

intron_exp_rmcentro <- subset(intron_exp, intron_exp$youngold!="Centro")

intron_exp_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_exp_compart_sep.txt', header = T)
str(intron_exp_sep)

intron_exp_sep_rmcentro <- subset(intron_exp_sep, intron_exp_sep$youngold!="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronmean_hypo_youngold.pdf", width=8, height=8)
ggplot(intron_exp_sep_rmcentro, aes(x=youngold, y=intron_mean_a1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dark grey","dodgerblue3"), labels=c("Long intron","Equal","Short intron"), name="Hypothecised direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,200) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Mean intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronnr_hypo_youngold.pdf", width=8, height=8)
ggplot(intron_exp_sep_rmcentro, aes(x=youngold, y=intron_nr_a1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dark grey","dodgerblue3"), labels=c("Long intron","Equal","Short intron"), name="Hypothecised direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()



#load the corresponding data files, on 21Jan2019
intron_total <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/21jan2019/Mvsl_a1a2_intron_exp_gencomp.txt', header = T)
str(intron_total)

intron_total_rmcentro <- subset(intron_total, intron_total$youngold!="Centro")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronmeanlength_diff_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=youngold, y=diffintronmean, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","dark grey","dodgerblue4")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-30,25) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_introntotallength_diff_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=youngold, y=diffintrontotal, fill=DE)) + 
  scale_fill_manual(values = c("firebrick4","dark grey","dodgerblue4")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-200,250) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronnr_diff_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=youngold, y=diffintronnr, fill=DE)) + 
  scale_fill_manual(values = c("firebrick4","dark grey","dodgerblue4")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-1,5)+                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#scatter plot 21jan.2019
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronmeanlength_diff_scatter_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=logFC.A1.A2, y=diffintronmean,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  ylim(-200,200) +
  labs(x='Log2(A1/A2)', y= 'Difference of average intron length (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_intronr_diff_scatter_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=logFC.A1.A2, y=diffintronnr,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  labs(x='Log2(A1/A2)', y= 'Difference of mean intron number (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_introntotallength_diff_scatter_youngold.pdf", width=8, height=8)
ggplot(intron_total_rmcentro, aes(x=logFC.A1.A2, y=diffintrontotal,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  labs(x='Log2(A1/A2)', y= 'Difference of total intron length (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm (logFC.A1.A2 ~ diffintronnr*DE-1, data=intron_total_rmcentro)
summary(y)
###
Estimate Std. Error t value Pr(>|t|)    
diffintronnr        0.562344   0.070340   7.995 1.47e-15 ***
  DEDown             -1.825360   0.031370 -58.188  < 2e-16 ***
  DENON              -0.007431   0.004710  -1.578    0.115    
DEUp                1.782527   0.034885  51.098  < 2e-16 ***
  diffintronnr:DENON -0.523949   0.080681  -6.494 8.83e-11 ***
  diffintronnr:DEUp  -0.476527   0.089739  -5.310 1.12e-07 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4219 on 8348 degrees of freedom
Multiple R-squared:  0.4238,	Adjusted R-squared:  0.4234 
F-statistic:  1023 on 6 and 8348 DF,  p-value: < 2.2e-16
###

y1 <- lm (logFC.A1.A2 ~ diffintronmean*DE-1, data=intron_total_rmcentro)
summary(y1)
###
Estimate Std. Error t value Pr(>|t|)    
diffintronmean       -0.006162   0.002709  -2.274  0.02298 *  
  DEDown               -1.837019   0.031554 -58.219  < 2e-16 ***
  DENON                -0.007488   0.004725  -1.585  0.11307    
DEUp                  1.774039   0.035232  50.353  < 2e-16 ***
  diffintronmean:DENON  0.006761   0.002803   2.412  0.01587 *  
  diffintronmean:DEUp   0.010276   0.003144   3.269  0.00108 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4233 on 8348 degrees of freedom
Multiple R-squared:   0.42,	Adjusted R-squared:  0.4196 
F-statistic:  1008 on 6 and 8348 DF,  p-value: < 2.2e-16
###

y2 <- lm (logFC.A1.A2 ~ diffintrontotal*DE-1, data=intron_total_rmcentro)
summary(y2)
###
Estimate Std. Error t value Pr(>|t|)    
diffintrontotal        0.0035823  0.0005771   6.207 5.65e-10 ***
  DEDown                -1.8143157  0.0315352 -57.533  < 2e-16 ***
  DENON                 -0.0074741  0.0047179  -1.584    0.113    
DEUp                   1.7832965  0.0351987  50.664  < 2e-16 ***
  diffintrontotal:DENON -0.0035044  0.0006426  -5.453 5.09e-08 ***
  diffintrontotal:DEUp  -0.0033309  0.0006430  -5.180 2.27e-07 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4226 on 8348 degrees of freedom
Multiple R-squared:  0.4219,	Adjusted R-squared:  0.4215 
F-statistic:  1015 on 6 and 8348 DF,  p-value: < 2.2e-16
###


###loading data for separate a1 and a2. 21jan.2019
intron_total_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/21jan2019/Mvsl_a1a2_intron_exp_gencomp_sep.txt', header = T)
str(intron_total_sep)

intron_total_sep_rmcentro <- subset(intron_total_sep,intron_total_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronnr_sep_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1intronnr, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronnr_sep_nobias_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1intronnr, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue4","firebrick4"), labels=c("A1","A2"), name="Haploid") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronmeanlength_sep_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1intronmean, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,250) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronmeanlength_sep_nobias_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1intronmean, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue4","firebrick4"), labels=c("A1","A2"), name="Haploid") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,220) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Avereage intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_introntotallength_sep_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1introntotal, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1500) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_introntotallength_sep_nobias_youngold.pdf", width=8, height=8)
ggplot(intron_total_sep_rmcentro, aes(x=youngold, y=a1introntotal, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue4","firebrick4"), labels=c("A1","A2"), name="Haploid") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1300) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Total intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

y4 <- lm(a1intronmean ~ haploid*youngold, data=intron_total_sep_rmcentro)
summary(y4) #intron mean length, only color strata makes siginificant difference. 

y5 <- lm(a1introntotal ~ haploid*youngold, data=intron_total_sep_rmcentro)
summary(y5) #no difference.

y6 <- lm(a1intronnr ~ haploid*youngold, data=intron_total_sep_rmcentro)
summary(y6) # no difference

#load the corresponding data files.
intron_total_fi <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/haploidwater_DEnonDE_intronfi08012019.txt', header = T)
str(intron_total_fi) #loading Jan.08.2019

intron_diff <- subset.data.frame(intron_total_fi, intron_total_fi$haploid == "A1")
str(intron_diff)

#plot intron number Jan.08.2019
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronnr_8genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total_fi, aes(x=chr, y=intronnr, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Avereage intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronnr_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total_fi, aes(x=chrom, y=intronnr, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Avereage intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronnrdiff_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_diff, aes(x=chrom, y=intronnrdiff, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.01,0.01) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Avereage intron number difference') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


###scatter plot on positive correlation between intron length and coding sequencing length.
ggplot(intron_total, aes(x=introntotal, y=genelength)) + geom_point()

#plot intron mean length
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronmeanlength_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total_fi, aes(x=chrom, y=meanlength, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
   geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,250) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Average intron length per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_intronmeanlengthdiff_8genomiccompartment.pdf", width=8, height=8)
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

ggplot(intron_diff, aes(x=chr, y=intronmeandiff, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  #stat_summary(fun.data = give.n, geom = "text",colour = "red", width = 1, size =5) +
  ylim(-10,10) +                    
  #scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Avereage intron length difference') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#plot intron total length
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_introntotallength_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total_fi, aes(x=chrom, y=totallength, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1","Not biased at A2","A1 biased at A1","A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1500) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Total intron length per gene') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_introntotallengthdiff_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_diff, aes(x=chrom, y=introntotaldiff, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Total intron length difference') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_meanintronlengthdiff_3genomiccompartment.pdf", width=8, height=8)
ggplot(intron_diff, aes(x=chrom, y=intronmeandiff, fill=interaction(haploid,DE_status))) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Average intron length difference') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
y1 <- lm(logFCA1A2~ intronmeandiff * intronnrdiff * introntotaldiff, data=intron_diff)
summary(y1)
##
lm(formula = logFCA1A2 ~ intronmeandiff * intronnrdiff * introntotaldiff, 
   data = intron_diff)

Residuals:
  Min       1Q   Median       3Q      Max 
-12.4884  -0.1419   0.0342   0.1888   7.3495 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)                                 -3.201e-02  5.874e-03  -5.449  5.2e-08 ***
  intronmeandiff                               5.308e-03  1.994e-03   2.662  0.00777 ** 
  intronnrdiff                                 8.277e-02  8.607e-02   0.962  0.33624    
introntotaldiff                              3.344e-04  8.305e-04   0.403  0.68721    
intronmeandiff:intronnrdiff                 -6.215e-03  1.908e-03  -3.257  0.00113 ** 
  intronmeandiff:introntotaldiff               1.524e-05  6.045e-06   2.521  0.01171 *  
  intronnrdiff:introntotaldiff                 1.119e-03  3.878e-04   2.886  0.00391 ** 
  intronmeandiff:intronnrdiff:introntotaldiff -2.494e-06  4.733e-06  -0.527  0.59827   
####

#### intron differences
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_singlecopy_intronnrdiffnormalized_genomiccompartment.pdf", width=8, height=8)
ggplot(intron_total, aes(x=Genomicloc, y=intronnrovergenelength)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.007) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron number difference between A1 and A2') +
  theme(axis.title.x = element_text(size=10,colour = "black"), axis.title.y = element_text(size=10,colour = "black")) +
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


##########################################################################################################################
######################  codes below are not really needed for final version, keep it just in case ########################
######################                                                                            ########################
###DE genes and intron size, number
intron_de <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/haploidwater_DEnonDE_intron_genomiccompart_a1a2sep_fi.txt', header = T)
str(intron_de)

ggplot(intron_de, aes(x=chr, y=Expression, fill=haploid)) + scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,10) +              
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Autosome", "PAR","Green", "Red","Orange","Black","Blue","Purple")) + 
  labs(x='Genomic compartment', y='Intron total length') +
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
