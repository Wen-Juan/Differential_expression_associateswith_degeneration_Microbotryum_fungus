#install R packages
library(gplots)
library(ggplot2)

#load datasets
datapath <- "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/"
outpath <- paste("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/", sub_analyse, sep="")


perc_count <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1_haploidwater_count_perct_matrix.txt",header = T)
str(perc_count)


ggplot(perc_count, aes(x=sin(perc), y=lqqogCPM,color=chr, shape=chr)) + 
  geom_point()+
  geom_smooth(method = lm)


ggplot(perc_count, aes(x=perc, y=logFC.A1.A2, color=chr, shape=chr)) + 
  geom_point() +
  geom_smooth()

##MAT
perc_count_MAT <- subset(perc_count,perc_count$chr == 'aMAT')
str(perc_count_MAT)

ggplot(perc_count_MAT, aes(x=sin(perc), y=logFC.A1.A2, color=chr, shape=chr)) + 
  geom_point() +
  geom_smooth(method=lm)

ggplot(perc_count_MAT, aes(x=sin(perc), y=lqqogCPM, color=chr, shape=chr)) + 
  geom_point() +
  geom_smooth(method=lm)

y <- lm(lqqogCPM ~ sin(perc)-1, perc_count)
anova(y)
summary(y)
