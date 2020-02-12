#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
install.packages("gridExtra")
library(gridExtra)
require(cowplot)
library(Hmisc)

#load the corresponding data files, between a1 and a2 homologs within Mvsl species, homolog length was trimmed to equal and removing TE genes. modified at jan.28.2019
dNdS_new <- read.table('~/input/dNdS/Mvsl_dnds_exp_gencomp.txt', header = T)
str(dNdS_new)

pdf("~/output/figures/Mvsl_dN_DEnonDe_twoparts28012019.pdf", width=8, height=8)
p_dN <- ggplot(dNdS_new, aes(x=youngold, y=dn, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_dN_updownnon_twoparts28012019.pdf", width=8, height=8)
p_dN3 <- ggplot(dNdS_new, aes(x=youngold, y=dn, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), guide = FALSE)  +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.7)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


pdf("~/output/figures/Mvsl_ds_DEnonDe_twoparts28012019.pdf", width=8, height=8)
p_dS <- ggplot(dNdS_new, aes(x=youngold, y=ds, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


pdf("~/output/figures/Mvsl_ds_updownnon_twoparts28012019.pdf", width=8, height=8)
p_dS3 <- ggplot(dNdS_new, aes(x=youngold, y=ds, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.08) +  
  theme_bw() + theme(legend.position = c(0.2, 0.7)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


#combine two figures in one
pdf("~/output/figures/Mvsl_dn_ds_updownnon_combine_new.pdf", width=14, height=8)
par(mar=c(5,5,4,3))
plot_grid(p_dN3, p_dS3, p_dnds3, labels=c('A','B','C'))
dev.off()

pdf("~/output/figures/Mvsl_dn_ds_combine_new.pdf", width=14, height=8)
par(mar=c(5,5,4,3))
plot_grid(p_dN, p_dS, labels=c('A','B'))
dev.off()

pdf("~/output/figures/Mvsl_dnds_DEnonDe_twoparts28012019.pdf", width=8, height=8)
 ggplot(dNdS_new, aes(x=youngold, y=dnds, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,7) +  
  theme_bw() + theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=14),axis.text.y = element_text(colour="black",size=14))
dev.off()

pdf("~/output/figures/Mvsl_dnds_DEnonDe3_twoparts28012019.pdf", width=8, height=8)
p_dnds3 <- ggplot(dNdS_new, aes(x=youngold, y=dnds, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("A2 biased","Not biased","A1 biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,7) +  
  theme_bw() + theme(legend.position = c(0.2, 0.7)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("~/output/figures/Mvsl_dn_ds_combine3_new.pdf", width=14, height=8)
par(mar=c(8,8,6,6))
plot_grid(p_dN3, p_dS3, p_dnds3, labels=c('A','B','C'))
dev.off()

dNdS_youngstrata <- subset(dNdS_new, dNdS_new$youngold=="ColorStrata")
mean(dNdS_youngstrata$dn[dNdS_youngstrata$DE2=='DE']) #0.013
mean(dNdS_youngstrata$dn[dNdS_youngstrata$DE2=='NON']) #0.005730769
mean(dNdS_youngstrata$ds[dNdS_youngstrata$DE2=='DE']) #0.0128

dNdS_oldstrata <- subset(dNdS_new, dNdS_new$youngold=="OldStrata")
mean(dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='DE']) #0.05336977
mean(dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='NON']) #0.03256509

dNdS_PAR <- subset(dNdS_new, dNdS_new$youngold=="bPAR")
mean(dNdS_PAR$dn[dNdS_PAR$DE2=='DE']) #0
mean(dNdS_PAR$dn[dNdS_PAR$DE2=='NON']) #0.000132

dNdS_auto <- subset(dNdS_new, dNdS_new$youngold=="Auto")
mean(dNdS_auto$dn[dNdS_auto$DE2=='DE']) #0
mean(dNdS_auto$dn[dNdS_auto$DE2=='NON']) #0


mean(dNdS_youngstrata$dnds[dNdS_youngstrata$DE2=='DE']) #1.787892
x<-dNdS_youngstrata$dnds[dNdS_youngstrata$DE2=='DE']
se <-function(x) sqrt(var(x)/length(x))
se(x)  #0.4606744
y <- na.omit(dNdS_youngstrata$dnds[dNdS_youngstrata$DE2=='NON'])
mean(y)  #1.477885
se <-function(x) sqrt(var(x)/length(x))
se(y)  
#0.3295856

y1 <- na.omit(dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='DE'])
mean(y1 ) # 2.369368
se <-function(x) sqrt(var(x)/length(x))
se(y1)  #0.192115
y2 <- na.omit(dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='NON'])
mean(y2) #2.297156
se <-function(x) sqrt(var(x)/length(x))
se(y2)  #0.1828743


wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$dn[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 3125, p-value = 0.0003954
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$ds[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 3136.5, p-value = 0.0003296
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='DE'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 2338, p-value = 0.3901

#stats, separating A1- and A2-biased genes.
dNdS_youngstrata <- subset(dNdS_new, dNdS_new$youngold=="ColorStrata")
dNdS_oldstrata <- subset(dNdS_new, dNdS_new$youngold=="OldStrata")

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1156.5, p-value = 0.01943
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1028.5, p-value = 0.1721

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1968.5, p-value = 0.002656
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 2108, p-value = 0.0001539

wilcox.test(dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dn[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 226, p-value = 0.8113
wilcox.test(dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$ds[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 154, p-value = 0.1223

wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1001, p-value = 0.05473
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Down'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='NON'],exact = FALSE) 
#W = 1337, p-value = 0.819
wilcox.test(dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Up'],dNdS_oldstrata$dnds[dNdS_oldstrata$DE=='Down'],exact = FALSE) 
#W = 281, p-value = 0.04061

###scatter point figures.
pdf("~/output/figures/Mvsl_dn_exp_correlation_DEpool_bw.pdf", width=8, height=8)
p_dncor <- ggplot(dNdS_new, aes(x=dn, y=abs,color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide=FALSE) + 
  geom_point(size=2.5) + geom_smooth(method = lm) + 
  theme_bw() + 
  theme(legend.position = c(0.4, 0.8)) +
  labs(x='dN', y='Gene expression ratio in |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_dn_ds_separatea1a2bias_new.pdf", width=14, height=8)
par(mar=c(8,8,6,6))
plot_grid(p_dncor1, p_dscor1, p_dndscor1, labels=c('A','B','C'))
dev.off()

pdf("~/output/figures/Mvsl_ds_exp_correlation_pool_bw.pdf", width=8, height=8)
p_dscor <- ggplot(dNdS_new, aes(x=ds, y=abs,color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) + 
  theme_bw() + xlim(0,0.25) +
  theme(legend.position = c(0.4, 0.8)) +
  labs(x='dS', y='Gene expression ratio in |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_dnds_exp_correlation_DEpool.pdf", width=14, height=8)
p_dnds <- ggplot(dNdS_new, aes(x=dnds, y=abs,color=DE2,shape=DE2)) +
  scale_shape_manual(values=c(16,1),labels=c("DE","Non-DE"), name = "") +
  scale_color_manual(values = c("black","dark grey"),guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8)) +
  xlim(0,12.5) +
  labs(x='dN/dS', y='Gene expression ratio in |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_dn_ds_combine_new_twofigures_bw.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(p_dncor, p_dscor, p_dnds,labels=c('A','B','C'))
dev.off()


pdf("~/output/figures/Mvsl_dnds_exp_correlation_DEnotpool.pdf", width=8, height=8)
ggplot(dNdS_new, aes(x=dnds, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(0,13) +
  labs(x='dN/dS', y='Absolute value of expression in Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

##load data of dN/dS between Mvsl and Mvsv
dNdS_MvslMvsv <- read.table('/input/MvSl_Mint/Mvsl_Mvsv_dnds_exp_all_fi.txt', header = T)
str(dNdS_MvslMvsv)


pdf("~/output/figures/MvslMvsv_dNdiffe_youngandold.pdf", width=8, height=8)
p1 <- ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dndiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.03,0.03) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvsv_dSdiffe_youngandold.pdf", width=8, height=8)
p2 <- ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dsdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),  guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.03,0.03) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvsv_dNdsdiffe_youngandold.pdf", width=8, height=8)
p3 <- ggplot(dNdS_MvslMvsv, aes(x=youngold, y=dndsdiff, fill=DE)) +
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("a2-bias","Not-bias","a1-bias"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-0.2,0.2) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Difference in dN/dS (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvld_dn_ds_combine_new_3figures_bw.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(p1, p2, p3,labels=c('A','B','C'))
dev.off()


dNdS_MvslMvsv_sep <- read.table('~/input/MvSl_Mint/Mvsl_Mvsv_dnds_exp_all_fi_sep.txt', header = T)
str(dNdS_MvslMvsv_sep)
dNdS_MvslMvsv_sep_oldstrata <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold=="OldStrata")
x <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Lowmutations")
y <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Vhighmutations")
z <- subset(dNdS_MvslMvsv_sep_oldstrata,dNdS_MvslMvsv_sep_oldstrata$DE2 == "Neutral")

wilcox.test(dNdS_MvslMvsv_sep$dna2[dNdS_MvslMvsv_sep$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep$dna2[dNdS_MvslMvsv_sep$DE2=='Lowmutations'],exact = FALSE) 
#W = 16941, p-value = 0.9902

mean(x$dna2) #0.03965833
mean(x$dsa2) #0.1263271
mean(x$dndsa2) # 0.3424813

mean(z$dna2) #0.03380762
mean(z$dsa2) #0.1370848
mean(z$dndsa2) #0.2839695

mean(y$dna2) #0.0394625
mean(y$dsa2) # 0.1274812
mean(y$dndsa2) # 0.336725

pdf("~/output/figures/MvslMvsv_dN_mutations_youngoldstrata.pdf", width=8, height=8)
pa <- ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dna2, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.45) + 
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvsv_ds_mutations_youngoldstrata.pdf", width=8, height=8)
pb <- ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dsa2, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,0.45) + 
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvsv_dNds_mutations_youngoldstrata.pdf", width=8, height=8)
pc <- ggplot(dNdS_MvslMvsv_sep, aes(x=youngold, y=dndsa2, fill=DE2)) + 
 scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"), labels=c("Low mutations","Equal mutations","High mutations"), name="Expectations") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,1) +  
  theme_bw() + 
  theme(legend.position = c(0.5, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata", "Old strata")) + 
  labs(x='Genomic compartment', y='dN/dS') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvld_dn_ds_combine_new_3figures_bw.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(pa, pb, pc,labels=c('A','B','C'))
dev.off()

dNdS_MvslMvsv_sep1 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "OldStrata")

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1267, p-value = 0.9094
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1262.5, p-value = 0.9341
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],exact = FALSE) 
#W = 1247.5, p-value = 0.989

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 4580.5, p-value = 0.1137
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5558, p-value = 0.6695
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Vhighmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 6059, p-value = 0.145

wilcox.test(dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 4474, p-value = 0.07169
wilcox.test(dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dna2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 5519, p-value = 0.7289
wilcox.test(dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Lowmutations'],dNdS_MvslMvsv_sep1$dndsa2[dNdS_MvslMvsv_sep1$DE2=='Neutral'],exact = FALSE) 
#W = 6045, p-value = 0.1531

dNdS_MvslMvsv_sep2 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "ColorStrata")

dNdS_MvslMvsv_sep3 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "Auto")
wilcox.test(dNdS_MvslMvsv_sep3$dsa2[dNdS_MvslMvsv_sep3$DE2=='Lowmutations'],dNdS_MvslMvsv_sep3$dsa2[dNdS_MvslMvsv_sep3$DE2=='Vhighmutations'],exact = FALSE) 
#W = 8060.5, p-value = 0.9952
wilcox.test(dNdS_MvslMvsv_sep3$dna2[dNdS_MvslMvsv_sep3$DE2=='Lowmutations'],dNdS_MvslMvsv_sep3$dna2[dNdS_MvslMvsv_sep3$DE2=='Vhighmutations'],exact = FALSE) 
#W = 8066, p-value = 0.9986
wilcox.test(dNdS_MvslMvsv_sep3$dndsa2[dNdS_MvslMvsv_sep3$DE2=='Lowmutations'],dNdS_MvslMvsv_sep3$dndsa2[dNdS_MvslMvsv_sep3$DE2=='Vhighmutations'],exact = FALSE) 
#W = 8057, p-value = 0.9905

dNdS_MvslMvsv_sep4 <- subset(dNdS_MvslMvsv_sep,dNdS_MvslMvsv_sep$youngold == "bPAR")
wilcox.test(dNdS_MvslMvsv_sep4$dsa2[dNdS_MvslMvsv_sep4$DE2=='Lowmutations'],dNdS_MvslMvsv_sep4$dsa2[dNdS_MvslMvsv_sep4$DE2=='Vhighmutations'],exact = FALSE) 
#W = 19, p-value = 0.936
wilcox.test(dNdS_MvslMvsv_sep4$dna2[dNdS_MvslMvsv_sep4$DE2=='Lowmutations'],dNdS_MvslMvsv_sep4$dna2[dNdS_MvslMvsv_sep4$DE2=='Vhighmutations'],exact = FALSE) 
#W = 18.5, p-value = 1
wilcox.test(dNdS_MvslMvsv_sep4$dndsa2[dNdS_MvslMvsv_sep4$DE2=='Lowmutations'],dNdS_MvslMvsv_sep4$dndsa2[dNdS_MvslMvsv_sep4$DE2=='Vhighmutations'],exact = FALSE) 
#W = 17, p-value = 0.936