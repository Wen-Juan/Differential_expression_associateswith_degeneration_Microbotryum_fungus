#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)

#load the corresponding data files. modify this on Jan.09.2019.
#de genes on MAT 
exp_DE_MAT_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_mat_cor.txt', header = T)
str(exp_DE_MAT_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_diffMAT_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_MAT_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.25,0.4) +                    
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
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.2,0.2) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#################################not consider this at the moment################################################
### richmedium
DE_richmedium <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_richmedium.txt', header = T)
str(DE_richmedium)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDEinterval_overlapwithTE_propgenes_richmeidum_cor.pdf", width=8, height=8)
ggplot(DE_richmedium, aes(x=interval, y=diffpropcor, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.1,0.1) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
#############################################################################################

### 1,20k up and down stream TE insertion sites among genomic compartments.
DE_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_2k10kupdownstream.txt', header = T)
str(DE_homolog)

################### I think this figure is not correct, both catogory makes the boxplot look difference among DE and Non-DE genes.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
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
ddev.off()
############################################################################################################
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff2kdownstream_logfc_3genomicparts", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-0.5,0.5) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff8kupstream_logfc_3genomicparts", width=8, height=8)
  ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
    labs(x='Genomic compartment') + 
    labs(y='Difference of TE number [2-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
  
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kdwonstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervaldowndiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
DE_homolog2 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_mod8.txt', header = T)
str(DE_homolog2)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20updiffa1,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20downdiffa2,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

