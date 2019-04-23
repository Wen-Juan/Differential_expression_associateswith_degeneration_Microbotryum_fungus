#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load dataset to randomdize the non-DE gene directions
nonDE_GC <- read.table('~/input/GCcontent/05feb2019/non-DE_GC.txt', header = T)
str(nonDE_GC)

GC0 <- data.frame(nonDE_GC$a1GC0, nonDE_GC$a2GC0)
GC3 <- data.frame(nonDE_GC$a1GC3, nonDE_GC$a2GC3)

GC0_rand <- randomizeMatrix(GC0,null.model = "richness",iterations = 1000)
GC3_rand <- randomizeMatrix(GC3,null.model = "richness",iterations = 1000)

total <- cbind(GC0_rand, GC3_rand)
head(total)

total_data <- cbind(nonDE_GC, total)
write.table(total_data,file = "~/input/GCcontent/Mvsl_a1a2_OGC_3GC_exp_genomic_low-high.txt",quote=F, row.names=T, sep='\t')

GC_ratio_rand <- read.table('~/input/GCcontent/Mvsl_a1a2_OGC_3GC_exp_genomic_low-high.txt', header = T)
str(GC_ratio_rand)

GC_ratio_rand_rmcentro <- subset(GC_ratio_rand, GC_ratio_rand$youngold != "Centro")

pdf("~/output/figures/Mvsl_a1a2_diffGC0_youngold.pdf", width=8, height=8)
paa1 <- ggplot(GC_ratio_rand_rmcentro, aes(x=youngold, y=-GC0diff, fill=DE2)) + 
  scale_fill_manual(values = c("white","dark grey"), guide = FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC0% difference') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

GC0_3_oldstrata <- subset(GC_ratio_rand_rmcentro, GC_ratio_rand_rmcentro$youngold == "OldStrata")
GC0_3_youngstrata <- subset(GC_ratio_rand_rmcentro, GC_ratio_rand_rmcentro$youngold == "ColorStrata")
GC0_3_auto <- subset(GC_ratio_rand_rmcentro, GC_ratio_rand_rmcentro$youngold == "Auto")
GC0_3_PAR <- subset(GC_ratio_rand_rmcentro, GC_ratio_rand_rmcentro$youngold == "bPAR")

wilcox.test(GC0_3_oldstrata$GC3diff[GC0_3_oldstrata$DE2=="DE"],GC0_3_oldstrata$GC3diff[GC0_3_oldstrata$DE2=="NON"], exact = FALSE)
#W = 4481, p-value = 0.8147
wilcox.test(GC0_3_oldstrata$GC0diff[GC0_3_oldstrata$DE2=="DE"],GC0_3_oldstrata$GC0diff[GC0_3_oldstrata$DE2=="NON"], exact = FALSE)
#W = 4136, p-value = 0.4988

wilcox.test(GC0_3_youngstrata$GC3diff[GC0_3_youngstrata$DE2=="DE"],GC0_3_youngstrata$GC3diff[GC0_3_youngstrata$DE2=="NON"], exact = FALSE)
#W = 67, p-value = 0.1737
wilcox.test(GC0_3_youngstrata$GC0diff[GC0_3_youngstrata$DE2=="DE"],GC0_3_youngstrata$GC0diff[GC0_3_youngstrata$DE2=="NON"], exact = FALSE)
#W = 47, p-value = 0.9244

wilcox.test(GC0_3_auto$GC3diff[GC0_3_auto$DE2=="DE"],GC0_3_auto$GC3diff[GC0_3_auto$DE2=="NON"], exact = FALSE)
#W = 1973300, p-value = 0.04557
wilcox.test(GC0_3_auto$GC0diff[GC0_3_auto$DE2=="DE"],GC0_3_auto$GC0diff[GC0_3_auto$DE2=="NON"], exact = FALSE)
#W = 1981300, p-value = 0.005744

wilcox.test(GC0_3_PAR$GC3diff[GC0_3_PAR$DE2=="DE"],GC0_3_PAR$GC3diff[GC0_3_PAR$DE2=="NON"], exact = FALSE)
#W = 594, p-value = 0.5601
wilcox.test(GC0_3_PAR$GC0diff[GC0_3_PAR$DE2=="DE"],GC0_3_PAR$GC0diff[GC0_3_PAR$DE2=="NON"], exact = FALSE)
#W = 566.5, p-value = 0.1924

pdf("~/output/figures/Mvsl_a1a2_diffGC3_youngold.pdf", width=8, height=8)
paa2 <- ggplot(GC_ratio_rand_rmcentro, aes(x=youngold, y=-GC3diff, fill=DE2)) + 
  scale_fill_manual(values = c("white","grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC3% difference') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/MvslMvld_GC_combine_new_2figures_bw.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(paa1, paa2,labels=c('A','B'))
dev.off()


pdf("~/output/figures/Mvsl_a1a2_diffGC3_corr_youngold.pdf", width=8, height=8)
paaa1 <- ggplot(GC_ratio_rand_rmcentro, aes(x=-GC3diff, y=abs,color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"),guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  xlim(-8,8) +
  ylim(0,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='GC3% difference', y= 'Gene expression ratio |Log2(A1/A2)|')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_a1a2_diffGC0_corr_youngold.pdf", width=8, height=8)
paaa2 <-ggplot(GC_ratio_rand_rmcentro, aes(x=-GC0diff, y=abs,color=DE2,shape=DE2)) +
  scale_shape_manual(values=c(16,1), labels=c("DE","Non-DE"), name = "Bias direction") +
  scale_color_manual(values = c("black","dark grey"), guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  xlim(-8,8) +
  ylim(0,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='GC0% difference', y= 'Gene expression ratio |Log2(A1/A2)|')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_GC_cor_combine_figs.pdf", width=8, height=6)
par(mar=c(8,8,6,6))
plot_grid(paaa1, paaa2, labels=c('A','B'))
dev.off()

cor.test(GC_ratio_rand_rmcentro$GC3diff[GC_ratio_rand_rmcentro$DE2 == "DE"], GC_ratio_rand_rmcentro$abs[GC_ratio_rand_rmcentro$DE2 == "DE"], method=c("pearson"))
#t = 1.1523, df = 592, p-value = 0.2496
cor.test(GC_ratio_rand_rmcentro$GC0diff[GC_ratio_rand_rmcentro$DE2 == "DE"], GC_ratio_rand_rmcentro$abs[GC_ratio_rand_rmcentro$DE2 == "DE"], method=c("pearson"))
#t = 1.6125, df = 592, p-value = 0.1074


###load data on 21Jan2019
GC_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/21jan2019/MvSl_a1a2_0GC_3GC_genomcomp_exp.txt', header = T)
str(GC_ratio)

GC_ratio_rmcentro <- subset(GC_ratio,GC_ratio$youngold !="Centro")
pdf("~/output/figures/Mvsl_a1a2_diffGC0_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=diffGC0, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_a1a2_diffGC3_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=diffGC3, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,4) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC3% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("~/output/figures/Mvsl_a1a2_GC3_scatter_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=logFC.A1.A2, y=diffGC3,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  labs(x='Log2(A1/A2)', y= 'Difference of %GC3 (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


##stats
wilcox.test(GC0_3$GC3diff[GC0_3$DE_status2=="Down"], GC0_3$GC3diff[GC0_3$DE_status2=="NON"], exact = FALSE)
#W = 2879, p-value = 0.3637
wilcox.test(GC0_3$GC3diff[GC0_3$DE_status2=="Up"], GC0_3$GC3diff[GC0_3$DE_status2=="NON"],
            exact = FALSE)
#W = 3377, p-value = 0.4419
wilcox.test(GC0_3$GC3diff[GC0_3$DE_status2=="Up"], GC0_3$GC3diff[GC0_3$DE_status2=="Down"], 
            exact = FALSE)
#W = 335.5, p-value = 0.2554