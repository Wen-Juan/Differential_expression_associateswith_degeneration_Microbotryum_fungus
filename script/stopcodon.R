#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
stopcodon <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/stopcodon/homolog_a1a2_cds_prot_length_removeTE_DEnonDEgeneexpree_intron_fi.txt', header = T)
str(stopcodon)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1A2_cdsratio_a1a2_3PARTS.pdf", width=8, height=8)
ggplot(stopcodon, aes(chr1,cdsratio)) + 
  #geom_boxplot() +
  geom_violin(fill = "blue") +
  #scale_x_discrete(labels=c("Autosome", "PAR","Green","Red", "Orange","Black","Blue","Purple")) + 
  scale_x_discrete(labels=c("Autosome", "PAR","NRR")) +
  labs(y='Coding sequence length ration (A1/A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_a1_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(stopcodon, aes(x=chr1, y=cdsratio, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey"), labels=c("Down","Up", "not bias"), name="DE expression") + 
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85) +
  #scale_x_discrete(labels=c("up:10-22k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()