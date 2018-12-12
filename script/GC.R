#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

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
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("Down","NON","Up"), name="Expression") + 
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  geom_jitter(aes(chr1_n + scat_adj, Godiff),
              position=position_jitter(width=0.05,height=0),
              alpha=0.4,
              size=1,
              show_guide=FALSE) +
  ylim(-2,2) +                    
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) + 
  labs(x='Genomic compartment', y='Difference in %GC between homologs (A1-A2)') +
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
  scale_fill_manual(values = c("firebrick2","grey","dodgerblue2"), labels=c("Down","NON","Up"), name="Expression") + 
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


