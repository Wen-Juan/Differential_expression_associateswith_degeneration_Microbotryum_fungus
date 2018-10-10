#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
compname <- 'tv46male'
contrast.name <- 'logFC.XY46.XX46'
group1 <- 'XY'
group2 <- 'XX'

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/LogCPM_0.05_tv46male_boxplot.txt",header = T)
#kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/LogCPM_0.05_tv46.txt", header = T)
#kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/LogCPM_0.05_tv43.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

#kdata <- subset(kdata,kdata$logFC.XY46.XX46<=-1)
#kdata <- subset(kdata,kdata$logFC.XY46.XX46>=1)

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- subset(map.data, map.data$chr=='Chr01') 
chr2.data <- subset(map.data, map.data$chr!='Chr01') 

group1.expr.minus <- (map.data$Tv1_463 + map.data$Tv2_461 + map.data$Tv2_464) / 3 #XY male group at G46
#group1.expr.minus <- (map.data$Tv1_432 + map.data$Tv3_431) / 2 #XY male group at G43

group2.expr.minus <- map.data$Tv2_462 #XX male group at G46
#group2.expr.minus <- (map.data$Tv5_431+map.data$Tv5_433+map.data$Tv6_432+map.data$Tv6_433)/4 #XX female group at G43
#group2.expr.minus <- (map.data$Tv1_462 + map.data$Tv5_461 + map.data$Tv6_461) / 3 #XX females at G46

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- log2(group1.expr/group2.expr)
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

#ggplot of boxplot
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_4a_g46.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-6,6)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Log2 ratio of XY0/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 

#########results of G43 stage
#W = 5180300, p-value = 0.08691
#########results of G46 stage
#W = 5116600, p-value = 0.1127

#calculate means by sliding windows
install.packages("RcppRoll")
library("RcppRoll")
install.packages("zoo")
library("zoo")
install.packages("microbenchmark")
library("microbenchmark")
require(zoo)


dn_ds <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds/amm88perc_annotation_dnds.txt", header = T)
str(dn_ds)

map.data <- dn_ds
chr1_1_start <-map.data$start[map.data$chr=='Chr02']
chr_1_end <-  180000000
  max(map.data$start[map.data$chr=='Chr02'])

zoo.dat <- zoo(map.data$dnds[map.data$chr=='Chr02'], c(chr1_1_start,chr_1_end))
y <- rollapply(zoo.dat, 40, FUN = mean, align = 'center', na.rm=TRUE) 

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/Figure_dnds_sexchr2.pdf", width=8, height=8)
par(mar=c(5,5,4,3)+0.6) 
#plot(c(0, chr_1_end), c(0,2), axes=F, lwd=2, xlab="Position (bp)", ylab="gene expression ratio log2(XY0/XX)", cex.axis=1.5, cex.lab=1.2, col="white")
plot(c(0, chr_1_end), c(0,0.4), axes=F, lwd=2, xlab="Position (bp)", ylab="dN/dS", cex.axis=1.5, cex.lab=1.2, col="white")
axis(1,c(0, 30000000,60000000,90000000, 120000000, 150000000, 180000000))
axis(2, c(0,0.2,0.4))
points(map.data$start[map.data$chr=='Chr02'], map.data$dnds[map.data$chr=='Chr02'], pch=20, lwd=2, type="p",col="gray50", main="",cex.axis=1.5)
#points(DC.data1dmrt$start[DC.data1dmrt$chr=='Chr01'],DC.data1dmrt$dN_dS[DC.data1dmrt$chr=='Chr01'], pch=20, lwd=3, type="p",col="red", main="",cex.axis=1.5)
lines(y,col="blue",lwd=4)
abline(v=116440117, col="blue",lwd=2,lty=3)
dev.off()

