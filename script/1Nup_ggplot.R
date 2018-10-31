
#install and load relavant R packages
install.packages("faraway")
install.packages("ggplot2")
library(faraway)
library(ggplot2)

#load the dataset
#compare the ratio of male biased gene among fold changes across stages
N_up <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A1allgenes_v1/LogCPM_logFC1_MAT_1Nup.txt", header = TRUE)
str(N_up)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC1_MAT_a1_1nupregulate.pdf", width=7,height=5)
ggplot(N_up, aes(x = start, y=logFC.water.diploid,group=share), cex=2) + 
  geom_point(size = 2, aes(color=share))+ 
  scale_color_manual(values=c("dark grey","red")) +
  labs (x = "Gene location on MAT chromosome", y= "Log2FC(N/N+N)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

A2N_up <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A2hapdi_v1/LogFC1_A2_1npregulate.txt", header = TRUE)
str(A2N_up)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC1_MAT_a2_1nupregulate.pdf", width=7,height=5)
ggplot(A2N_up, aes(x = start, y=logFC.water.diploid, group=share), cex=2) + 
  geom_point(size = 2, aes(color=share))+ 
  scale_color_manual(values=c("dark grey","red")) +
  labs (x = "Gene location on MAT chromosome", y= "Log2FC(N/N+N)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

#A1 log2FC intron mean length

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_intron_mean_exp_cor.pdf", width=7,height=5)
A1_LogFC_de_nonde <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_DE__nonDE_intron_cor.txt", header = TRUE)
str(A1_LogFC_de_nonde)

Abovezero_LogFC <- subset(A1_LogFC_de_nonde, log2FC >= 0)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_up_NONde_intron_mean_exp_cor.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(mean), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_mean)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
  #+
  #geom_smooth(method='lm',aes(fill=group))

dev.off()

model3 <- lm(abs(log2FC) ~ mean*group, A1_LogFC_de_nonde)
anova(model3) #both intron number and DE has significant inflence.

###
Response: abs(log2FC)
Df Sum Sq Mean Sq   F value Pr(>F)    
mean          1  46.06   46.06  394.6083 <2e-16 ***
  group         1 635.58  635.58 5444.7625 <2e-16 ***
  mean:group    1   0.27    0.27    2.3538  0.125    
Residuals  8065 941.44    0.12                     
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
###

#intron total length
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_all_intron_sum_expression_cor.pdf", width=7,height=5)

ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(sum), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_sum)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) +
  geom_smooth()

dev.off()
model4 <- lm(abs(log2FC) ~ sum*group, A1_LogFC_de_nonde)
 #both intron number and DE has significant inflence.
anova(model4)
###
Response: abs(log2FC)
Df Sum Sq Mean Sq  F value    Pr(>F)    
sum          1 230.27  230.27 1987.919 < 2.2e-16 ***
  group        1 452.47  452.47 3906.081 < 2.2e-16 ***
  sum:group    1   6.39    6.39   55.134 1.241e-13 ***
  Residuals 8065 934.23    0.12  

###

###intron size
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_intron.pdf", width=7,height=5)
A1_LogFC_de_nonde <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_DE__nonDE_intron_cor.txt", header = TRUE)
str(A1_LogFC_de_nonde)

Abovezero_LogFC <- subset(A1_LogFC_de_nonde, log2FC >= 0)

ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(intron_nr), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_number)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) +
  geom_smooth(method='lm',aes(fill=group))

dev.off()

model1 <- lm(log2FC ~ intron_nr*group, A1_LogFC_de_nonde)
anova(model1) #both intron number and DE has significant inflence.

###
Response: log2FC
Df  Sum Sq Mean Sq F value    Pr(>F)    
intron_nr          1    6.23  6.2294  23.094 1.570e-06 ***
  group              1   15.99 15.9851  59.262 1.544e-14 ***
  intron_nr:group    1   21.53 21.5278  79.810 < 2.2e-16 ***
  Residuals       8065 2175.44  0.2697    


### to incorporate hemizygous genes, need to use the normalized expression level:
A1_LogFC_de_nonde_count <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_logFC_count_DE_nonDE_hemi_intron_cor.txt", header = TRUE)
str(A1_LogFC_de_nonde_count)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(y = ave_count, x=log2(mean), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (x = "Log2(intron_mean)", y= "Normalized counts", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

model6 <- lm(ave_count ~ mean*group, A1_LogFC_de_nonde_count)
anova(model6)

##
Response: ave_count
Df Sum Sq Mean Sq F value    Pr(>F)    
mean          1    113  113.27 24.2203 8.758e-07 ***
  group         2    919  459.68 98.2897 < 2.2e-16 ***
  mean:group    2      3    1.67  0.3576    0.6994    
Residuals  8266  38659    4.68 
##

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_number.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(y = ave_count, x=log2(intron_nr), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.4)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (x = "Log2(intron_number)", y= "Normalized counts", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

model7 <- lm(ave_count ~ log2(intron_nr)*group, A1_LogFC_de_nonde_count)
anova(model7)

##
Response: ave_count
Df Sum Sq Mean Sq  F value  Pr(>F)    
log2(intron_nr)          1   3692  3692.2 866.1445 < 2e-16 ***
  group                    2    744   371.8  87.2232 < 2e-16 ***
  log2(intron_nr):group    2     22    11.2   2.6183 0.07298 .  
Residuals             8266  35236     4.3  
##

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_sum.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(y = ave_count, x=log2(sum), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.4)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (x = "Log2(intron_sum)", y= "Normalized counts", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

model8 <- lm(ave_count ~ log2(sum)*group, A1_LogFC_de_nonde_count)
anova(model8)

##
Response: ave_count
Df Sum Sq Mean Sq  F value    Pr(>F)    
log2(sum)          1      8    7.52   1.6246    0.2025    
group              2   1294  647.07 139.7748 < 2.2e-16 ***
  log2(sum):group    2    127   63.32  13.6773 1.174e-06 ***
  Residuals       8266  38266    4.63 
##
