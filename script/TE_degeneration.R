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
