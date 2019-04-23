#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load data to reandomdize the non-DE gene expression direction
nonDE_intron <- read.table('~/input/intron_degeneration/05feb2019/intron_random_diff.txt', header = T)
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
write.table(total_data,file = "~/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_compart_randomdized.txt",quote=F, row.names=T, sep='\t')

intron_random <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/intron_degeneration/05feb2019/Mvsl_a1a2_intron_exp_compart_low-high.txt', header = T)
str(intron_random)

intron_random_rmcentro <- subset(intron_random, intron_random$youngold != "Centro")

pdf("~/output/figures/Mvsl_introntotal_diff1_youngold.pdf", width=8, height=8)
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


pdf("~/output/figures/Mvsl_intronmean_diff1_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_intronnumber_diff1_youngold.pdf", width=8, height=8)
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


pdf("~/output/figures/Mvsl_intron_combine_3figs_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(pa,pb,pc,labels=c('A','B','C'))
dev.off()


pdf("~/output/figures/Mvsl_a1a2_intronnr_corr_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_a1a2_intronmeanlength_corr_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_a1a2_introntotallength_corr_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_intron_3figures.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(pa1, pa3, pa2,labels=c('A','B','C'))
dev.off()


####check coding sequence length and protein lengh variation across genomic compartment
cds_prot <- read.table('~/input/intron_degeneration/Mvsl_a1a2_cds_prot_compart.txt', header = T)
str(cds_prot)

cds_prot_rmcentro <- subset(cds_prot, cds_prot$youngold != "Centro")

pdf("~/output/figures/Mvsl_cds_combined_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_prot_combined_youngold.pdf", width=8, height=8)
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

pdf("~/output/figures/Mvsl_coding_proteinratio_combine_figs.pdf", width=12, height=8)
par(mar=c(8,8,6,6))
plot_grid(pa1, pa2, labels=c('A','B'))
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


intron_exp <- read.table('~/input/intron_degeneration/Mvsl_a1a2_intron_exp_compart.txt', header = T)
str(intron_exp)

intron_exp_rmcentro <- subset(intron_exp, intron_exp$youngold!="Centro")

intron_exp_sep <- read.table('~/input/intron_degeneration/Mvsl_a1a2_intron_exp_compart_sep.txt', header = T)
str(intron_exp_sep)

intron_exp_sep_rmcentro <- subset(intron_exp_sep, intron_exp_sep$youngold!="Centro")

pdf("~/output/figures/Mvsl_intronmean_hypo_youngold.pdf", width=8, height=8)
ggplot(intron_exp_sep_rmcentro, aes(x=youngold, y=intron_mean_a1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dark grey","dodgerblue3"), labels=c("Long intron","Equal","Short intron"), name="Hypothecised direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,200) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Mean intron length') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("~/output/figures/Mvsl_intronnr_hypo_youngold.pdf", width=8, height=8)
ggplot(intron_exp_sep_rmcentro, aes(x=youngold, y=intron_nr_a1, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","dark grey","dodgerblue3"), labels=c("Long intron","Equal","Short intron"), name="Hypothecised direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,15) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Intron number') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()