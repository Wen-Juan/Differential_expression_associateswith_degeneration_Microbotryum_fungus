#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load dataset to randomdize the non-DE gene directions, on 15feb.2019.
nonDE_GC <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/05feb2019/non-DE_GC.txt', header = T)
str(nonDE_GC)

GC0 <- data.frame(nonDE_GC$a1GC0, nonDE_GC$a2GC0)
GC3 <- data.frame(nonDE_GC$a1GC3, nonDE_GC$a2GC3)

GC0_rand <- randomizeMatrix(GC0,null.model = "richness",iterations = 1000)
GC3_rand <- randomizeMatrix(GC3,null.model = "richness",iterations = 1000)

total <- cbind(GC0_rand, GC3_rand)
head(total)

total_data <- cbind(nonDE_GC, total)
write.table(total_data,file = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/05feb2019/Mvsl_a1a2_GC_nonDE_compart_randomdized.txt",quote=F, row.names=T, sep='\t')

GC_ratio_rand <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/05feb2019/Mvsl_a1a2_OGC_3GC_exp_genomic_rand.txt', header = T)
str(GC_ratio_rand)

GC_ratio_rand_rmcentro <- subset(GC_ratio_rand, GC_ratio_rand$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_youngold.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC3_youngold.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvld_GC_combine_new_2figures_bw.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(paa1, paa2,labels=c('A','B'))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC3_corr_youngold.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_corr_youngold.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_GC_cor_combine_figs.pdf", width=8, height=6)
par(mar=c(8,8,6,6))
plot_grid(paaa1, paaa2, labels=c('A','B'))
dev.off()



y1 <- lm(abs ~ DE2/GC3diff-1, data=GC_ratio_rand_rmcentro)
summary (y1)

######
Estimate Std. Error t value Pr(>|t|)    
DE2DE          1.890582   0.018109 104.402  < 2e-16 ***
  DE2NON         0.214062   0.004947  43.267  < 2e-16 ***
  DE2DE:GC3diff  0.099984   0.025561   3.912 9.24e-05 ***
  DE2NON:GC3diff 0.001061   0.020745   0.051    0.959 
Multiple R-squared:  0.5992,	Adjusted R-squared:  0.599 
F-statistic:  3193 on 4 and 8544 DF,  p-value: < 2.2e-16
#######

y2 <- lm(abs ~ DE2/GC0diff-1, data=GC_ratio_rand_rmcentro)
summary (y2)

#########
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE          1.890261   0.018088 104.504  < 2e-16 ***
  DE2NON         0.214033   0.004942  43.305  < 2e-16 ***
  DE2DE:GC0diff  0.338532   0.061858   5.473 4.56e-08 ***
  DE2NON:GC0diff 0.027840   0.041860   0.665    0.506 
Multiple R-squared:  0.5999,	Adjusted R-squared:  0.5997 
F-statistic:  3203 on 4 and 8544 DF,  p-value: < 2.2e-16
##########




#load the corresponding data files, on 05Feb.2019
GC_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/05feb2019/Mvsl_a1a2_OGC_3GC_exp_genomic.txt', header = T)
str(GC_ratio)

GC_ratio_rmcentro <- subset(GC_ratio, GC_ratio$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=GC0diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC0% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC3_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=GC3diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC3% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_corr_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=GC0diff, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-5,5) +
  ylim(0,13) +
  labs(x='Difference of %GC0 (A1-A2)', y= 'Absolute value of ratio Log2(A1/A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC3_corr_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=GC3diff, y=abs,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-6,6) +
  ylim(0,13) +
  labs(x='Difference of %GC3 (A1-A2)', y= 'Absolute value of ratio Log2(A1/A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(abs~GC3diff*DE-1, data=GC_ratio_rmcentro)
y1 <- lm(abs~GC3diff+DE, data=GC_ratio_rmcentro)
anova(y,y1)
summary(y)
summary(y1)

y2 <- lm(abs~GC0diff*DE-1, data=GC_ratio_rmcentro)
y3 <- lm(abs~GC0diff+DE-1, data=GC_ratio_rmcentro)
anova(y2,y3)
summary(y2)
summary(y3)

GC_ratio_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/05feb2019/Mvsl_a1a2_OGC_3GC_exp_genomic_sep.txt', header = T)
str(GC_ratio_sep)

GC_ratio_sep_rmcentro <- subset(GC_ratio_sep, GC_ratio_sep$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0GC3_youngold_compat.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC0, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"),labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1", "Not biased at A2","A1 biased at A1", "A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(45,70) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC0 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC3_youngold_compat.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC3, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"),labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1", "Not biased at A2","A1 biased at A1", "A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(40,90) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC3 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0_expexted_youngold_compat.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC0, fill=DE2)) + 
  scale_fill_manual(values = c("dodgerblue3" ,"grey", "firebrick3"),labels=c("High GC","Equal", "Low GC"), name="Expected GC direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(45,70) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC0 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC3_expected_youngold_compat.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC3, fill=DE2)) + 
  scale_fill_manual(values = c("dodgerblue3" ,"grey", "firebrick3"),labels=c("High GC","Equal", "Low GC"), name="Expected GC direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(40,90) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC3 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()



#load the corresponding data files, on 21Jan.2019
GC_ratio_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/21jan2019/MvSl_a1a2_0GC_3GC_genomcomp_exp_sep.txt', header = T)
str(GC_ratio_sep)

GC_ratio_sep_rmcentro <- subset(GC_ratio_sep,GC_ratio_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC0, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"),labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1", "Not biased at A2","A1 biased at A1", "A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(45,70) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC0 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC3_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_sep_rmcentro, aes(x=youngold, y=a1GC3, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"),labels=c("A2 biased at A1","A2 biased at A2","Not biased at A1", "Not biased at A2","A1 biased at A1", "A1 biased at A2"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(40,85) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Percentage of GC3 on coding sequence') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###load data on 21Jan2019
GC_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/21jan2019/MvSl_a1a2_0GC_3GC_genomcomp_exp.txt', header = T)
str(GC_ratio)

GC_ratio_rmcentro <- subset(GC_ratio,GC_ratio$youngold !="Centro")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=diffGC0, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2.5,2.5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC3_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=youngold, y=diffGC3, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,4) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='GC3% difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC3_scatter_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=logFC.A1.A2, y=diffGC3,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  labs(x='Log2(A1/A2)', y= 'Difference of %GC3 (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()
y <- lm(logFC.A1.A2~diffGC3*DE-1, data=GC_ratio_rmcentro)
summary(y)
############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
diffGC3       -0.100181   0.032919  -3.043 0.002348 ** 
  DEDown        -1.835780   0.031489 -58.299  < 2e-16 ***
  DENON         -0.007271   0.004725  -1.539 0.123914    
DEUp           1.787692   0.034810  51.356  < 2e-16 ***
  diffGC3:DENON  0.138499   0.036991   3.744 0.000182 ***
  diffGC3:DEUp   0.084836   0.043565   1.947 0.051526 .  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4232 on 8348 degrees of freedom
Multiple R-squared:  0.4202,	Adjusted R-squared:  0.4198 
F-statistic:  1008 on 6 and 8348 DF,  p-value: < 2.2e-16
###############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_diffGC0_scatter_youngold.pdf", width=8, height=8)
ggplot(GC_ratio_rmcentro, aes(x=logFC.A1.A2, y=diffGC0,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  xlim(-13,13) +
  labs(x='Log2(A1/A2)', y= 'Difference of %GC3 (A1-A2)')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y1 <- lm(logFC.A1.A2~diffGC0*DE-1, data=GC_ratio_rmcentro)
summary(y1)
###
lm(formula = logFC.A1.A2 ~ diffGC0 * DE - 1, data = GC_ratio_rmcentro)

Residuals:
  Min       1Q   Median       3Q      Max 
-10.7066  -0.1445   0.0203   0.1724   5.5261 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
diffGC0        0.023937   0.067182   0.356  0.72162    
DEDown        -1.831093   0.031487 -58.154  < 2e-16 ***
  DENON         -0.007314   0.004725  -1.548  0.12167    
DEUp           1.791334   0.034801  51.473  < 2e-16 ***
  diffGC0:DENON  0.050901   0.075929   0.670  0.50264    
diffGC0:DEUp  -0.302251   0.109208  -2.768  0.00566 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4232 on 8348 degrees of freedom
Multiple R-squared:  0.4202,	Adjusted R-squared:  0.4198 
F-statistic:  1008 on 6 and 8348 DF,  p-value: < 2.2e-16
####

GC_ratio_mat <- subset(GC_ratio_rmcentro, GC_ratio_rmcentro$youngold !="Auto")
str(GC_ratio_mat)

y2 <- lm(logFC.A1.A2~diffGC0*DE-1, data=GC_ratio_mat)
summary(y2)
#####
lm(formula = logFC.A1.A2 ~ diffGC0 * DE - 1, data = GC_ratio_mat)

Residuals:
  Min      1Q  Median      3Q     Max 
-9.8354 -0.2796 -0.0163  0.2677  4.2945 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
diffGC0       -0.05431    0.16359  -0.332    0.740    
DEDown        -2.64591    0.17944 -14.745   <2e-16 ***
  DENON          0.05639    0.06323   0.892    0.373    
DEUp           2.14658    0.15398  13.941   <2e-16 ***
  diffGC0:DENON  0.17881    0.19714   0.907    0.365    
diffGC0:DEUp  -0.34936    0.29104  -1.200    0.231    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.025 on 335 degrees of freedom
Multiple R-squared:  0.5535,	Adjusted R-squared:  0.5455 
F-statistic: 69.21 on 6 and 335 DF,  p-value: < 2.2e-16
#####

y3 <- lm(logFC.A1.A2~diffGC3*DE-1, data=GC_ratio_mat)
summary(y3)
#####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
diffGC3       -0.10959    0.08533  -1.284   0.1999    
DEDown        -2.64985    0.17897 -14.806   <2e-16 ***
  DENON          0.06076    0.06353   0.956   0.3396    
DEUp           2.11720    0.15321  13.819   <2e-16 ***
  diffGC3:DENON  0.17366    0.10332   1.681   0.0937 .  
diffGC3:DEUp   0.13603    0.11311   1.203   0.2300    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.027 on 335 degrees of freedom
Multiple R-squared:  0.5519,	Adjusted R-squared:  0.5439 
F-statistic: 68.76 on 6 and 335 DF,  p-value: < 2.2e-16
#####

###loading data on GC and dnds on 21Jan.2019
GC_dnds <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/21jan2019/Mvsl_exp_GC_dnds.txt', header = T)
str(GC_dnds)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Cor_a1GC3_dn_scatter.pdf", width=8, height=8)
ggplot(GC_dnds, aes(x=a1GC3, y=dn, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,0.25) +
 # xlim(-13,13) +
  labs(x='%GC3of A1', y= 'dN')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Cor_a1GC0_dn_scatter.pdf", width=8, height=8)
ggplot(GC_dnds, aes(x=a1GC0, y=dn, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,0.25) +
  # xlim(-13,13) +
  labs(x='%GC0 of A1', y= 'dN')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Cor_a2GC3_dn_scatter.pdf", width=8, height=8)
ggplot(GC_dnds, aes(x=a2GC3, y=dn, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,0.25) +
  # xlim(-13,13) +
  labs(x='%GC3 of A2', y= 'dN')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Cor_a2GC0_dn_scatter.pdf", width=8, height=8)
ggplot(GC_dnds, aes(x=a2GC0, y=dn, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,0.25) +
  # xlim(-13,13) +
  labs(x='%GC0 of A2', y= 'dN')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Cor_a2GC0_dn_scatter.pdf", width=8, height=8)
ggplot(GC_dnds, aes(x=ds, y=diffGC0, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  #ylim(0,0.25) +
  # xlim(-13,13) +
  #labs(x='%GC0 of A1', y= 'dS')+
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load the corresponding data files.
GC0 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/Mvsl_a1a2_DEnonDE_intron_GCcodinggenes.txt', header = T)
str(GC0)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0_8compartments_withdots.pdf", width=8, height=8)
GC0$chr1_n[GC0$chr1 == "aAutosome"] <- 1
GC0$chr1_n[GC0$chr1 == "bPAR"] <- 2
GC0$chr1_n[GC0$chr1 == "NRR"] <- 3

GC0$scat_adj[GC0$DE_status2 == "Down"] <- -0.25
GC0$scat_adj[GC0$DE_status2 == "NON"] <- 0
GC0$scat_adj[GC0$DE_status2 == "Up"] <- 0.25
ggplot(GC0, aes(x=chr1, y=Godiff, fill=DE_status2)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  #geom_jitter(aes(chr1_n + scat_adj, Godiff),
             # position=position_jitter(width=0.05,height=0),
             # alpha=0.4,
              #size=1,
             # show_guide=FALSE) +
  ylim(-2,2) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference in total %GC between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###stats
GC_NRR <- subset(GC0, GC0$chr1 == "NRR")

wilcox.test(GC_NRR$Godiff[GC_NRR$DE_status2=="Down"], GC_NRR$Godiff[GC_NRR$DE_status2=="NON"], alternative = "greater",
            exact = FALSE, correct = FALSE)
#W = 2545.5, p-value = 0.5584
wilcox.test(GC_NRR$Godiff[GC_NRR$DE_status2=="Up"], GC_NRR$Godiff[GC_NRR$DE_status2=="NON"], alternative = "greater",
            exact = FALSE, correct = FALSE)
#W = 3694.5, p-value = 0.4768
wilcox.test(GC_NRR$Godiff[GC_NRR$DE_status2=="Up"], GC_NRR$Godiff[GC_NRR$DE_status2=="Down"], alternative = "greater",
            exact = FALSE, correct = FALSE)
#W = 418.5, p-value = 0.4342

######
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0_3compartments.pdf", width=8, height=8)
ggplot(GC0, aes(x=chr1, y=Godiff, fill=DE_status2)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-1,1) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference in %GC between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1_GC0_3compartments.pdf", width=8, height=8)

GC0$chr1_n[GC0$chr1 == "aAutosome"] <- 1
GC0$chr1_n[GC0$chr1 == "bPAR"] <- 2
GC0$chr1_n[GC0$chr1 == "NRR"] <- 3

GC0$scat_adj[GC0$DE_status2 == "Down"] <- -0.25
GC0$scat_adj[GC0$DE_status2 == "NON"] <- 0
GC0$scat_adj[GC0$DE_status2 == "Up"] <- 0.25

ggplot(GC0, aes(x=chr1, y=GC0A1, fill=DE_status2)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("Down","NON","Up"), name="Expression") + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,alpha=0.85) +
  geom_jitter(aes(chr1_n + scat_adj, GC0A2),
              position=position_jitter(width=0.05,height=0),
              alpha=0.4,
              size=1,
              show_guide=FALSE) +
  ylim(45,65) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='%GC in A1') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###GC0 and GC3
GC0_3 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/GCcontent/Mvsl_a1a2_DEnonDE_intron_GC0_3_codinggenes.txt', header = T)
str(GC0_3)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2_GC0_3_3compartments_withdots.pdf", width=8, height=8)
GC0_3$chr1_n[GC0_3$chr1 == "aAutosome"] <- 1
GC0_3$chr1_n[GC0_3$chr1 == "bPAR"] <- 2
GC0_3$chr1_n[GC0_3$chr1 == "NRR"] <- 3

GC0_3$scat_adj[GC0_3$DE_status2 == "Down"] <- -0.25
GC0_3$scat_adj[GC0_3$DE_status2 == "NON"] <- 0
GC0_3$scat_adj[GC0_3$DE_status2 == "Up"] <- 0.25
ggplot(GC0_3, aes(x=chr1, y= GC3diff, fill=DE_status2)) + 
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("Down","NON","Up"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  geom_jitter(aes(chr1_n + scat_adj, Godiff),
              position=position_jitter(width=0.05,height=0),
              alpha=0.4,
              size=1,
              show_guide=FALSE) +
  ylim(-2,2) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference in %GC3 between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
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