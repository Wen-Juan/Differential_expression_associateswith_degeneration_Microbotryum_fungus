
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

ggplot(A2N_up, aes(x = start, y=logFC.water.diploid, group=share), cex=2) + 
  geom_boxplot(size = 2, aes(color=share))+ 
  scale_color_manual(values=c("dark grey","red")) +
  labs (x = "Gene location on MAT chromosome", y= "Log2FC(N/N+N)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

#A1 log2FC intron mean length
A1_LogFC_de_nonde <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_DE__nonDE_intron_cor_sort_fi.txt", header = TRUE)
str(A1_LogFC_de_nonde)

Abovezero_LogFC <- subset(A1_LogFC_de_nonde, log2FC >= 0)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_up_NONde_intron_mean_exp_cor_corrected.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(mean), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_mean)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
  #geom_smooth(method='lm',aes(fill=group))

dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_up_NONde_intron_mean_exp_cor_corrected_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde, aes(x=log2FC, y=log2(mean), fill=group)) + scale_fill_manual(values = c("firebrick2","grey40"),labels=c("DE genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,varwidth = TRUE) +
  ylim(5.2,8) +
  labs(y="Log2(Intron_mean_length)") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
    theme(legend.position = "right")

dev.off()

model3 <- lm(abs(log2FC) ~ mean*group, A1_LogFC_de_nonde)
anova(model3) #both intron number and DE has significant inflence.

###
Response: abs(log2FC)
Df Sum Sq Mean Sq   F value Pr(>F)    
mean          1   8.50    8.50   72.9547 <2e-16 ***
  group         1 675.44  675.44 5798.9282 <2e-16 ***
  mean:group    1   0.02    0.02    0.2136  0.644    
Residuals  8065 939.39    0.12 
###

#intron total length
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_all_intron_sum_expression_cor_corrected.pdf", width=7,height=5)

ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(sum), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6,varwidth = TRUE)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_sum)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
#+
 # geom_smooth()
dev.off()


pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_all_intron_sum_expression_cor_corrected_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde, aes(x=log2FC, y=sum, fill=group)) + scale_fill_manual(values = c("firebrick2","grey40"),labels=c("DE genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,varwidth = TRUE) +
  ylim(200,1300) +
  labs(y="Intron_total_length") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")

ggplot(A1_LogFC_de_nonde, aes(x=log2FC, y=log2(sum), fill=group)) +
  scale_fill_manual(values = c("firebrick2","grey40"),labels=c("DE genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,varwidth = TRUE) +
  ylim(6.3,11) +
  labs(y="Log2(Intron_total_length)") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")

dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_all_intron_sum_expression_cor_corrected_boxplot.pdf", width=7,height=5)

ggplot(A1_LogFC_de_nonde, aes(x=log2FC, y=log2(sum), fill=group)) + scale_fill_manual(values = c("firebrick2","grey40"),labels=c("DE genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,varwidth = TRUE) +
  ylim(7,10) +
  labs(y="Intron_total_length") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
dev.off()


model4 <- lm(abs(log2FC) ~ log2(sum)*group, A1_LogFC_de_nonde)
 #both intron number and DE has significant inflence.
anova(model4)
###
Response: abs(log2FC)
Df Sum Sq Mean Sq  F value    Pr(>F)    
log2(sum)          1   3.68    3.68   31.771 1.793e-08 ***
  group              1 681.25  681.25 5873.676 < 2.2e-16 ***
  log2(sum):group    1   3.02    3.02   26.060 3.386e-07 ***
  Residuals       8065 935.40    0.1
###

###intron size
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_intron_corrected.pdf", width=7,height=5)
A1_LogFC_de_nonde <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_DE__nonDE_intron_cor_sort_fi.txt", header = TRUE)
str(A1_LogFC_de_nonde)

Abovezero_LogFC <- subset(A1_LogFC_de_nonde, log2FC >= 0)

ggplot(A1_LogFC_de_nonde, aes(y = log2FC, x=log2(intron_nr), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (x = "Log2(intron_number)", y= "Log2FC(A1/A2)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
#+
 # geom_smooth(method='lm',aes(fill=group))
dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_intron_corrected_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde, aes(x=log2FC, y=intron_nr, fill=group)) + scale_fill_manual(values = c("firebrick2","grey40"),labels=c("DE genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,varwidth = TRUE) +
  ylim(0,15) +
  labs(y="Intron number per gene") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
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
A1_LogFC_de_nonde_count <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_a1_logFC_count_DE_nonDE_hemi_intron_cor_sort_fi.txt", header = TRUE)
str(A1_LogFC_de_nonde_count)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(y = log2(mean), x=ave_count, group=group), cex=2)  + 
  geom_point(size = 2, aes(color=group),alpha = 0.8)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (x = "Log2(intron_mean)", y= "Normalized counts", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_mean_cor_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x=ave_count, y=log2(mean), fill=group)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"),labels=c("DE genes","Hemizygous genes","non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,position = position_dodge2(preserve = "single")) +
  ylim(5.5,8.5) +
  labs(y="Log2(Intron_mean_length)") +
  scale_x_discrete(labels=c("DE genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
dev.off()

model6 <- lm(ave_count ~ mean*group, A1_LogFC_de_nonde_count)
anova(model6)

##
Analysis of Variance Table

Response: ave_count
Df Sum Sq Mean Sq  F value Pr(>F)    
mean          1      1    0.88   0.1891 0.6636    
group         2   1032  515.81 110.2884 <2e-16 ***
  mean:group    2      3    1.30   0.2789 0.7566    
Residuals  8266  38660    4.68  
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

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_number_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x=ave_count, y=intron_nr, fill=group)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"),labels=c("DE genes","Hemizygous genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,position = position_dodge2(preserve = "single")) +
  ylim(1,10) +
  labs(y="Intron number per gene") +
  scale_x_discrete(labels=c("DE genes","Hemizygous genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
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
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_sum_boxplot.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x=ave_count, y=log(sum), fill=group)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"),labels=c("DE genes","Hemizygous genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,position = position_dodge2(preserve = "single")) +
  ylim(4,8) +
  labs(y="Log2(Intron total length per gene)") +
  scale_x_discrete(labels=c("DE genes","Hemizygous genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC0_DE_NONde_hemi_normalizedcount_intron_sum_boxplot2.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x=ave_count, y=sum, fill=group)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"),labels=c("DE genes","Hemizygous genes","Non-DE genes")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.4,alpha=0.8,position = position_dodge2(preserve = "single")) +
  ylim(0,1100) +
  labs(y="Log2(Intron total length per gene)") +
  scale_x_discrete(labels=c("DE genes","Hemizygous genes","Non-DE genes"),name="Log2FC(A1/A2)") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  theme(legend.position = "right")
dev.off()

model8 <- lm(ave_count ~ log2(sum)*group, A1_LogFC_de_nonde_count)
anova(model8)

##
Response: ave_count
Df Sum Sq Mean Sq  F value Pr(>F)    
log2(sum)          1   3454  3453.9 806.8903 <2e-16 ***
  group              2    849   424.7  99.2095 <2e-16 ***
  log2(sum):group    2      9     4.6   1.0817 0.3391    
Residuals       8266  35382     4.3  
##

## correlation of intron number and mean size
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/DE_NONde_hemi_normalizedcount_A1intron_sum_vd_mean.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x = log2(mean), y=log2(sum)), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (y = "Log2(intron_total_length)", x= "Log2(intron_mean)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/DE_NONde_hemi_normalizedcount_A1intron_number_vs_mean.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x = intron_nr, y=log2(mean), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (y = "Log2(intron_total length)", x= "Intron number", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/DE_NONde_hemi_normalizedcount_A1intron_number_vs_totallength.pdf", width=7,height=5)
ggplot(A1_LogFC_de_nonde_count, aes(x = intron_nr, y=log2(sum), group=group), cex=2) + 
  geom_point(size = 2, aes(color=group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","blue","grey")) +
  labs (y = "Log2(intron_total length)", x= "Intron number", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()


###
a1a2_LogFC_de_nonde_count <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/water_haploid_de_nonde_a1_a2_list_intron_fi.txt", header = TRUE)
str(a1a2_LogFC_de_nonde_count)

ggplot(a1a2_LogFC_de_nonde_count, aes(x = ave_count, y=log2(intron_a1sum), group=Group), cex=2) + 
  geom_point(size = 2, aes(color=Group),alpha = 0.6)+ 
  scale_color_manual(values=c("red","grey")) +
  labs (y = "Log2(intron_total length)", x= "Intron number", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()
