#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load data on 13feb2019 for directional ratio, and randomdize non-DE genes directions for ratio calculation.
diff_indel <- read.table('~/input/stopcodon/04Feb2019/prot_diff_genomic_compartment.txt', header = T)
str(diff_indel)

pdf("~/output/figures/Mvsl_early_stopcodon_protdiff.pdf", width=8, height=8)
pa <- ggplot(diff_indel, aes(x=Comp, y=prop, fill=DE)) + 
  scale_fill_manual(values = c("white","grey"),labels=c("DE","Non-DE"), name="Bias direction") +
  geom_bar(position="dodge", stat="identity", alpha=0.85, width=0.6, color = "black") +
  ylim(0,0.7) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome+PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion of genes with different stop codon positions') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


##stats
old_de <- c(25,18)
old_nonde <- c(56,75)
prop.test(old_de, old_nonde) #X-squared = 5.2951, df = 1, p-value = 0.02138

auto_de <- c(4,7)
auto_nonde <- c(40,21)
prop.test(auto_de, auto_nonde) #X-squared = 3.6165, df = 1, p-value = 0.05721


old_de <- c(25,7)
auto_de <- c(56,21)
prop.test(old_de,auto_de) #X-squared = 0.40605, df = 1, p-value = 0.524

old_nonde <- c(18,4)
auto_nonde <- c(75,40)
prop.test(old_nonde,auto_nonde) #X-squared = 2.462, df = 1, p-value = 0.1166

diff_indel2 <- read.table('~/input/stopcodon/diff_indel.txt', header = T)
str(diff_indel2)

pdf("~/output/figures/Mvsl_indels_mean_compartment.pdf", width=8, height=8)
pb <- ggplot(diff_indel2, aes(x=comp, y=indel, fill=DE)) +
  scale_fill_manual(values = c("white","grey"),guide =FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85, color = "black") +
  ylim(0,7) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome+PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Mean indel number') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_earlystopcodon_indels_combine2figs.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(pb,pa,labels=c('A','B'))
dev.off()

diff_indel2_oldstrata <- subset(diff_indel2, diff_indel2$comp =="old")
diff_indel2_auto <- subset(diff_indel2, diff_indel2$comp =="Auto")
wilcox.test(diff_indel2_oldstrata$indel[diff_indel2_oldstrata$DE=="DE"], diff_indel2_oldstrata$indel[diff_indel2_oldstrata$DE=="NON"]) #W = 2453.5, p-value = 0.01331
wilcox.test(diff_indel2_auto$indel[diff_indel2_auto$DE=="DE"], diff_indel2_auto$indel[diff_indel2_auto$DE=="NON"]) #W = 490.5, p-value = 0.02548

non_DE_protein <- read.table('~/input/stopcodon/Mvsl_a1a2_proteinlength_nonDE.txt', header = T)
str(non_DE_protein)

prot_ratio <- data.frame(non_DE_protein$a1prot, non_DE_protein$a2prot)
cds_ratio <- data.frame(non_DE_protein$a1cds, non_DE_protein$a2cds)

prot_ratio_rand <- randomizeMatrix(prot_ratio,null.model = "richness",iterations = 1000)
cds_ratio_rand <- randomizeMatrix(cds_ratio,null.model = "richness",iterations = 1000)

total_data <- cbind(non_DE_protein, prot_ratio_rand,cds_ratio_rand)
write.table(total_data,file = "~/input/stopcodon/Mvsl_a1a2_exp_gencompt_protlength_fi_higher-low.txt",quote=F, row.names=T, sep='\t')

diff_prot_length_rand <- read.table('~/input/stopcodon/04Feb2019/Mvsl_a1a2_exp_gencompt_protlength_fi_higher-low.txt', header = T)
str(diff_prot_length_rand1)

pdf("~/output/figures/Mvsl_protenlength_ratio_rand_youngold.pdf", width=8, height=8)
ggplot(diff_prot_length_rand1, aes(x=ratioprot, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(1,16),labels=c("Non-DE","DE"), name = "") +
  scale_color_manual(values = c("dark grey","black"),guide=FALSE) +
  geom_point(size=2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(0,2) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='Ratio of protein length between alleles', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

diff_prot_length_rand_rmcentro <- subset(diff_prot_length_rand,diff_prot_length_rand$youngold!="Centro") 
diff_prot_length_rand_oldstrata <- subset(diff_prot_length_rand,diff_prot_length_rand$youngold=="OldStrata")
diff_prot_length_rand_youngstrata <- subset(diff_prot_length_rand,diff_prot_length_rand$youngold=="ColorStrata")
diff_prot_length_rand_PAR <- subset(diff_prot_length_rand,diff_prot_length_rand$youngold=="bPAR")
diff_prot_length_rand_autosome <- subset(diff_prot_length_rand,diff_prot_length_rand$youngold=="Auto")
  

pdf("~/output/figures/Mvsl_protein_ratio_youngold_poolDE.pdf", width=8, height=8)
ggplot(diff_prot_length_rand_rmcentro, aes(x=youngold, y=ratioprot, fill=DE2)) +
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.75,1.25) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Protein length ratio') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("~/output/figures/Mvsl_cds_ratio_youngold_poolDE.pdf", width=8, height=8)
ggplot(diff_prot_length_rand_rmcentro, aes(x=youngold, y=ratiocds, fill=DE2)) +
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0.75,1.25) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Coding sequence length ratio') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

wilcox.test(diff_prot_length_rand_oldstrata$ratiocds[diff_prot_length_rand_oldstrata$DE2=='DE'],diff_prot_length_rand_oldstrata$ratiocds[diff_prot_length_rand_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4336, p-value = 0.8807
wilcox.test(diff_prot_length_rand_oldstrata$ratioprot[diff_prot_length_rand_oldstrata$DE2=='DE'],diff_prot_length_rand_oldstrata$ratioprot[diff_prot_length_rand_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4457, p-value = 0.8616


pdf("~/output/figures/Mvsl_cdslength_ratio_rand_youngold.pdf", width=8, height=8)
ggplot(diff_prot_length_rand, aes(x=ratiocds, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(0,2) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  labs(x='Ratio of coding sequence length', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


#load data on 12feb2019
diff_prot_length <- read.table('~/input/stopcodon/diff_protein_length.txt', header = T)
str(diff_prot_length)

pdf("~/output/figures/Mvsl_diff_protenlength_ratio_youngold.pdf", width=8, height=8)
ggplot(data=diff_prot_length, aes(x=genomcom,y=prop,fill=factor(DE))) +
  scale_fill_manual(values = c("white","dark grey"), labels=c("DE","Non-DE"), name="Bias") + 
  ylim(0,0.8) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.8)) +
  geom_bar(position="dodge",stat="identity",width=0.5,alpha=0.8,colour="black") + 
  geom_text(aes(label=total),position=position_dodge(width=0.6), vjust=-0.2) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(y='Proportion of genes with different protein length', x='Genomic compartments') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


####move orange strata genes to old strata.
stopcodon_ratio_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/19jan2019/Mvsl_a1a2_exp_cds_protein_compart_sep.txt', header = T)
str(stopcodon_ratio_sep)

stopcodon_ratio_sep1 <- subset(stopcodon_ratio_sep, stopcodon_ratio_sep$youngold != "Centro")
str(stopcodon_ratio_sep1)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_codingsequencea1a2_devidprotein3times_8PARTS.pdf", width=8, height=8)
ggplot(stopcodon_ratio_sep1, aes(x=youngold, y=cdsa2expest, fill=interaction(haploid,DE))) +
         scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2")) +
         geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
         ylim(0.975,1.05) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Ratio of coding sequence length/(3*protein length)') +
  theme_bw() + 
  theme(legend.position = c(0.85, 0.75)) +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###proportion of genes with length difference.
gene_ratio <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/19jan2019/number_gene_ratio.txt', header = T)
str(gene_ratio)

gene_ratio_cds <- subset(gene_ratio, gene_ratio$category =="cds")
gene_ratio_prot <- subset(gene_ratio, gene_ratio$category =="prot")

pdf("~/output/figures/Mvsl_ratiogene_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_cds, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","dodgerblue4","light grey","dark grey"), labels=c("DE & ratio!=1","DE & ratio=1","Non-DE & ratio!=1","Non-DE & ratio=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiogeneprot_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_prot, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","dodgerblue4","light grey","dark grey"), labels=c("DE & ratio!=1","DE & ratio=1","Non-DE & ratio!=1","Non-DE & ratio=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

gene_ratio_cds_notequal <- subset(gene_ratio_cds, gene_ratio_cds$type =="cdsnot")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_ratiocsd_youngold.pdf", width=8, height=8)
ggplot(gene_ratio_cds_notequal, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE & ratio!=1","Non-DE & ratio!=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()

gene_ratio_prot_notequal <- subset(gene_ratio_prot, gene_ratio_prot$type =="proteinanot")
pdf("~/output/figures/Mvsl_ratioprot_youngold1.pdf", width=8, height=8)
ggplot(gene_ratio_prot_notequal, aes(factor(Comp), ratio, fill = interaction(type, bias))) + 
  scale_fill_manual(values = c("dodgerblue2","light grey"), labels=c("DE & ratio!=1","Non-DE & ratio!=1"), name="Expression") + 
  geom_bar(stat="identity", position = "dodge",lpha=0.9,lwd=0.5) +
  geom_text(
    aes(label = value), 
    position = position_dodge(0.9),
    vjust = -0.4, size = 3.5
  ) +
  scale_x_discrete(labels=c("Autosome", "PAR","Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Proportion')
dev.off()
