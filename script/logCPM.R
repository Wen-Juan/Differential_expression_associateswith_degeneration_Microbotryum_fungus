#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load data of LogCPM  
logCPM_sig <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_differentialexpression_Microbotryumfungi/input/LogCPM/de_0.05_0_haploidwater_new4_copy.txt', header = T)
str(logCPM_sig)

logCPM_sig_a1 <- subset(logCPM_sig, logCPM_sig$logFC.A1.A2 > 0)
logCPM_sig_a1

logCPM_sig_a2 <- subset(logCPM_sig, logCPM_sig$logFC.A1.A2 < 0)
logCPM_sig_a2


ggplot(logCPM_sig, aes(chr, logCPM, fill = factor(bias)))+
  geom_boxplot()

logCPM_all <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_differentialexpression_Microbotryumfungi/input/LogCPM/LogCPM_0.05_haploidwater_new4.txt', header = T)
str(logCPM_all)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum_original_22042019_notongithub/output/figures/Mvsl_logCPM_all_chr.pdf", width=10, height=8)
ggplot(logCPM_all, aes(chr, logCPM, fill = factor(bias)))+
  scale_fill_manual(values = c("grey", "firebrick3","dodgerblue3"),labels=c("not bias","a1 bias","a2 bias"))  +
  geom_boxplot(notch=FALSE,alpha=0.8, color = "black") +
  ylim(0,15.5) +  
  theme_bw() + 
  theme(legend.position = c(0.85, 0.85)) +
  scale_x_discrete(labels=c("MAT","Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Others")) + 
  labs(x='Chromosomes', y='LogCPM') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

logCPM_all_short <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_differentialexpression_Microbotryumfungi/input/LogCPM/LogCPM_0.05_haploidwater_new4_short2.txt', header = T)
str(logCPM_all_short)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_differentialexpression_Microbotryumfungi/output/figures/Mvsl_logTPM_all_chr1.pdf", width=10, height=8)
ggplot(logCPM_all_short, aes(x=chr, y=logTPM, fill=interaction(matingtype,bias))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("a1-bias at a1","a1-bias at a2","No-bias at a1","No-bias at a2","a2-bias at a1","a2-bias at a2")) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,20) +  
  theme_bw() + 
  theme(legend.position = c(0.85, 0.8)) +
  scale_x_discrete(labels=c("MAT","Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Others")) + 
  labs(x='Chromosomes', y='Log2TPM') +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()
