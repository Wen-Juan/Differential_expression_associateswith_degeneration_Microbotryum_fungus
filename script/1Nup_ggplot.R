
#install and load relavant R packages
install.packages("faraway")
install.packages("ggplot2")
library(faraway)
library(ggplot2)

#load the dataset
#compare the ratio of male biased gene among fold changes across stages
N_up <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A1allgenes_v1/LogCPM_logFC1_MAT_1Nup.txt", header = TRUE)
str(N_up)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC1_MAT_a1_1nupregulate.pdf", width=7,height=5)
ggplot(N_up, aes(x = start, y=logFC.water.diploid,group=share), cex=2) + 
  geom_point(size = 2, aes(color=share))+ 
  scale_color_manual(values=c("dark grey","red")) +
  labs (x = "Gene location on MAT chromosome", y= "Log2FC(N/N+N)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()

A2N_up <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A2hapdi_v1/LogFC1_A2_1npregulate.txt", header = TRUE)
str(A2N_up)

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Log2FC1_MAT_a2_1nupregulate.pdf", width=7,height=5)
ggplot(A2N_up, aes(x = start, y=logFC.water.diploid, group=share), cex=2) + 
  geom_point(size = 2, aes(color=share))+ 
  scale_color_manual(values=c("dark grey","red")) +
  labs (x = "Gene location on MAT chromosome", y= "Log2FC(N/N+N)", size =2) +
  theme(axis.text.x = element_text(size = 12,color = "black"))  + 
  theme(axis.text.y = element_text(size = 12,color = "black")) +
  theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) 
dev.off()
