#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
te_interval <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/interval.txt', header = T)
str(te_interval)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_interval_overlapwithTE_genes.pdf", width=8, height=8)
ggplot(te_interval, aes(x=interval, y=numberofgenes, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_bar(stat="identity",position=position_dodge()) +
  ylim(0,1900) +                    
  scale_x_discrete(labels=c("up:0-2k", "up:2-10k","up:10-20k","down:0-2k", "down:2-10k","down:10-20k")) + 
  labs(y='Number of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2_interval_overlapwithTE_prop_genes.pdf", width=8, height=8)
ggplot(te_interval, aes(x=interval, y=prop, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue2","grey"), labels=c("A1","A2"), name="Haploid type") + 
  geom_bar(stat="identity",position=position_dodge()) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:0-2k", "up:2-10k","up:10-20k","down:0-2k", "down:2-10k","down:10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#differential expression and TE insertion site
###this is wrong
#exp_DE <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_wrong.txt', header = T)
#str(exp_DE)

#pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1_DEnonDEinterval_overlapwithTE_prop_genes.pdf", width=8, height=8)
#ggplot(exp_DE, aes(x=interval, y=prop, fill=bias)) + 
#  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
#  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
#  ylim(0,0.5) +                    
#  scale_x_discrete(labels=c("up:0-2k", "2-10k","10-20k","down:0-2k", "2-10k","10-20k")) + 
#  labs(y='Proportion of genes with TE insertion site') +
#  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
#  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
#dev.off()

#pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A1A2diff_DEnonDEinterval_overlapwithTE_prop_genes.pdf", width=8, height=8)
#ggplot(exp_DE, aes(x=interval, y=diffprop, fill=bias)) + 
  #scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  #geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  #ylim(-0.2,0.4) +                    
  #scale_x_discrete(labels=c("up:0-2k", "2-10k","10-20k","down:0-2k", "2-10k","10-20k")) + 
 # labs(y='Difference in prop. (A1-A2) of genes with TE insertion site') +
  #theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  #theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
#dev.off()

##A2 genome
exp_DE_a2 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_a2.txt', header = T)
str(exp_DE_a2)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_A2_DEnonDEinterval_overlapwithTE_prop_genes.pdf", width=8, height=8)
ggplot(exp_DE_a2, aes(x=interval, y=prop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(0,0.3) +                    
  scale_x_discrete(labels=c("up:0-2k", "2-10k","10-20k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

##this is correct
exp_DE_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_cor.txt', header = T)
str(exp_DE_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_cor, aes(x=interval, y=prop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(0,0.2) +                    
  scale_x_discrete(labels=c("up:10-22k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
  dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1a2diff_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_cor, aes(x=interval, y=diffpropcor, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(-0.03,0.05) +                    
  scale_x_discrete(labels=c("up:10-22k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#de genes on MAT 
exp_DE_MAT_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_mat_cor.txt', header = T)
str(exp_DE_MAT_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_diffMAT_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_MAT_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(-0.15,0.2) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

dev.off()

#DE autosome
exp_DE_auto_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_auto_cor.txt', header = T)
str(exp_DE_auto_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_autosomes_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_auto_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(-0.1,0.1) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### richmedium
DE_richmedium <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_richmedium.txt', header = T)
str(DE_richmedium)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDEinterval_overlapwithTE_propgenes_richmeidum_cor.pdf", width=8, height=8)
ggplot(DE_richmedium, aes(x=interval, y=diffpropcor, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("A1 bias","A2 bias", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  ylim(-0.02,0.02) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### 1:1 homolog
DE_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_2kupdownstream.txt', header = T)
str(DE_homolog)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_8genomicparts_2kupstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog, aes(y=logFC.A1.A2, x=k2updiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~chr) +
  ylim(-4,4) +  
  xlim(-2,2) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kdownstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog, aes(y=logFC.A1.A2, x=k2downdiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~genomiccomp) +
  ylim(-4,4) +  
  xlim(-2,2) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()

##1:1 homolog 2k10k upstream and downstream TE insertion sites.
DE_homolog_2k10k <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_2k10kupdownstream.txt', header = T)
str(DE_homolog_2k10k)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2k10kupstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog_2k10k, aes(y=logFC.A1.A2, x=k10kupstreamdiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~genomiccomp) +
  ylim(-2,2) +  
  xlim(-4,4) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2k10kdownstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog_2k10k, aes(y=logFC.A1.A2, x=k10kdownstreamdiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~genomiccomp) +
  ylim(-2,2) +  
  xlim(-4,4) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kupstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog_2k10k, aes(y=logFC.A1.A2, x=k8kintervalupdiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~genomiccomp) +
  ylim(-2,2) +  
  xlim(-4,4) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kdownstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog_2k10k, aes(y=logFC.A1.A2, x=k8kintervaldowndiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~genomiccomp) +
  ylim(-2,2) +  
  xlim(-4,4) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
dev.off()
