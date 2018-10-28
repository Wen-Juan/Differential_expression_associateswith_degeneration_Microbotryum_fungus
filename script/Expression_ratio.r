#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.
col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A2hapdi/'
kdata <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A2hapdi/LogCPM_0.05_A2hapdi.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

map.data <- subset(kdata, kdata$start!='NA')
#chr2.data <- subset(map.data, map.data$chr!='aMAT') 

group1.expr.minus <- (map.data$water1 + map.data$water2 ) / 2
group2.expr.minus <- (map.data$di1 + map.data$di2) / 2

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- log2(group1.expr/group2.expr)
#map.data$ratio <- group1.expr/group2.expr
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA
chr1.data <- data.frame(subset(map.data, map.data$chr=='aPR' | map.data$chr=='bHD1' | map.data$chr=='cHD2')) 
chr2.data <- data.frame(subset(map.data, map.data$chr!='aPR')) 
chr3.data <- data.frame(subset(chr2.data, chr2.data$chr!='bHD1')) 
chr4.data <- data.frame(subset(chr3.data, chr3.data$chr!='cHD2')) 
#
write.csv (map.data, "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A1MintMvSl/Logcpm_A1_mintmvsl.csv")

#gene expression ratio Log2(XY/XX)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MintMvsl_log2ratio.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","red","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-0.6,0.6)) + 
  geom_boxplot(notch = FALSE,outlier.shape=NA,)+
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of A1water/dikaryon gene expression') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#stats
##
wilcox.test(chr1.data$ratio,chr4.data$ratio,exact = FALSE) 
##W = 1596900, p-value = 1.542e-05

#calculate means by sliding windows
install.packages("RcppRoll")
library("RcppRoll")
install.packages("zoo")
library("zoo")
install.packages("microbenchmark")
library("microbenchmark")
require(zoo)
install.packages("rtracklayer")
library("rtracklayer")

map.data <- subset(map.data, map.data$ratio!='NaN')
##Assigning Transcripts to Chromosomes###
pr <- subset(map.data, map.data$chr=="aPR")
hd1 <- subset(map.data, map.data$chr=="bHD1")
hd2 <- subset(map.data, map.data$chr=="cHD2")
Chrmat <- subset(map.data, map.data$chr=="aMAT")
Chr1 <- subset(map.data, map.data$chr=="Chr01")
Chr2 <- subset(map.data, map.data$chr=="Chr02")
Chr3 <- subset(map.data, map.data$chr=="Chr03")
Chr4 <- subset(map.data, map.data$chr=="Chr04")
Chr5 <- subset(map.data, map.data$chr=="Chr05")
Chr6 <- subset(map.data, map.data$chr=="Chr06")
Chr7 <- subset(map.data, map.data$chr=="Chr07")
Chr8 <- subset(map.data, map.data$chr=="Chr08")
Chr9 <- subset(map.data, map.data$chr=="Chr09")
Chr10 <- subset(map.data, map.data$chr=="Chr10")
Chr11 <- subset(map.data, map.data$chr=="Chr11")
Chr12 <- subset(map.data, map.data$chr=="Chr12")
Chr13 <- subset(map.data, map.data$chr=="Chr13")
Chr14 <- subset(map.data, map.data$chr=="Chr14")
Chr15 <- subset(map.data, map.data$chr=="Chr15")
Chr16 <- subset(map.data, map.data$chr=="Chr16")
Chr17 <- subset(map.data, map.data$chr=="Chr17")
#Chr18 <- subset(map.data, map.data$chr=="Chr18")

## Sort According to Position ####
pr_sort <- pr[order(pr$start),]
hd1_sort <- hd1[order(hd1$start),]
hd2_sort <- hd2[order(hd2$start),]
Chrmat_sort <- Chrmat[order(Chrmat$start),] 
Chr01_sort <- Chr1[order(Chr1$start),] 
Chr02_sort <- Chr2[order(Chr2$start),] 
Chr03_sort <- Chr3[order(Chr3$start),] 
Chr04_sort <- Chr4[order(Chr4$start),] 
Chr05_sort <- Chr5[order(Chr5$start),] 
Chr06_sort <- Chr6[order(Chr6$start),]
Chr07_sort <- Chr6[order(Chr7$start),]
Chr08_sort <- Chr6[order(Chr8$start),]
Chr09_sort <- Chr6[order(Chr9$start),]
Chr10_sort <- Chr6[order(Chr10$start),]
Chr11_sort <- Chr6[order(Chr11$start),]
Chr12_sort <- Chr6[order(Chr12$start),]
Chr13_sort <- Chr6[order(Chr13$start),]
Chr14_sort <- Chr6[order(Chr14$start),]
Chr15_sort <- Chr6[order(Chr15$start),]
Chr16_sort <- Chr6[order(Chr16$start),]
Chr17_sort <- Chr6[order(Chr17$start),]
#Chr18_sort <- Chr6[order(Chr18$start),]


###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

prRT<- rollmean(smooth(na.approx(pr_sort$ratio)),10)
hd1RT<- rollmean(smooth(na.approx(hd1_sort$ratio)),10)
hd2RT<- rollmean(smooth(na.approx(hd2_sort$ratio)),10)
ChrmatRT<- rollmean(smooth(na.approx(Chrmat_sort$ratio)),10)
Chr1RT<- rollmean(smooth(na.approx(Chr01_sort$ratio)),10)
Chr2RT<- rollmean(smooth(na.approx(Chr02_sort$ratio)),10)
Chr3RT<- rollmean(smooth(na.approx(Chr03_sort$ratio)),10)
Chr4RT<- rollmean(smooth(na.approx(Chr04_sort$ratio)),10)
Chr5RT<- rollmean(smooth(na.approx(Chr05_sort$ratio)),10)
Chr6RT<- rollmean(smooth(na.approx(Chr06_sort$ratio)),10)
Chr6RT<- rollmean(smooth(na.approx(Chr08_sort$ratio)),10)
Chr7RT<- rollmean(smooth(na.approx(Chr08_sort$ratio)),20)
Chr8RT<- rollmean(smooth(na.approx(Chr09_sort$ratio)),10)
Chr9RT<- rollmean(smooth(na.approx(Chr11_sort$ratio)),10)
Chr10RT<- rollmean(smooth(na.approx(Chr11_sort$ratio)),20)
Chr11RT<- rollmean(smooth(na.approx(Chr12_sort$ratio)),10)
Chr12RT<- rollmean(smooth(na.approx(Chr14_sort$ratio)),10)
Chr13RT<- rollmean(smooth(na.approx(Chr14_sort$ratio)),20)
Chr14RT<- rollmean(smooth(na.approx(Chr16_sort$ratio)),10)
Chr15RT<- rollmean(smooth(na.approx(Chr16_sort$ratio)),20)
Chr16RT<- rollmean(smooth(na.approx(Chr17_sort$ratio)),10)
Chr17RT<- rollmean(smooth(na.approx(Chr18_sort$ratio)),10)

Rt<- c(ChrmatRT, Chr1RT,Chr2RT,Chr3RT,Chr4RT,Chr5RT,Chr6RT,Chr8RT,Chr9RT,Chr11RT,Chr12RT,Chr14RT,Chr16RT,Chr16RT)

myfunction <- function(i){
  Info <- sample(i,1,replace=FALSE)
  return(Info)
}

my.perm <- c()
for(i in 1:10^3){ my.perm[i] <- myfunction(Rt) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]

RMpalette <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")

#mating type chromosome
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/1N_NN_A2_log2_mat.pdf", width=7,height=5)

Chr_pos <-  rollmean(smooth(Chrmat_sort$start),10)
Chr_ratio<- rollmean(smooth(na.approx(Chrmat_sort$ratio)),10)
data.frame(Chr_pos,Chr_ratio)
plot(Chrmat_sort$start, Chrmat_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(N/N+N)",main="MAT genes")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2,col="red")
abline(h=highCI,lty=2,col="red")
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")

dev.off()

#mating type chromosome HD1
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/1N_NN_a2_Chr05_log2.pdf", width=7,height=5)

Chr_pos <-  rollmean(smooth(Chr05_sort$start),10)
Chr_ratio<- rollmean(smooth(na.approx(Chr05_sort$ratio)),10)
data.frame(Chr_pos,Chr_ratio)
plot(Chr05_sort$start, Chr05_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(N/N+N)",main="Chr05 genes")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2,col="red")
abline(h=highCI,lty=2,col="red")


dev.off()

#mating type chromosome hd2
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvSlMint_hd2_5genes_log2line.pdf", width=7,height=5)

Chr_pos <-  rollmean(smooth(hd2_sort$start),5)
Chr_ratio<- rollmean(smooth(hd2_sort$ratio),5)
data.frame(Chr_pos,Chr_ratio)
plot(hd2_sort$start, hd2_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(Mvsl/Mint)",main="HD2 genes")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")

dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/1N_NN_Chr05_log2.pdf", width=7,height=5)

Chr_pos <-  rollmean(smooth(Chr05_sort$start),10)
Chr_ratio<- rollmean(smooth(na.approx(Chr05_sort$ratio)),10)
data.frame(Chr_pos,Chr_ratio)
plot(Chr05_sort$start, Chr05_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(N/N+N)",main="Chr05 genes")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2,col="red")
abline(h=highCI,lty=2,col="red")
dev.off()
#autosome chromosome1
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMint_Chr01_sliding20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr01_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr01_sort$ratio)),20)
plot(Chr01_sort$start, Chr01_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-5,5), xlab="Position(bp)", ylab="Log2(MvSl/Mint)",main="Chr01")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()

##autosome chromosome2
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMint_Chr02_ratio_20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr02_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr02_sort$ratio)),20)
plot(Chr02_sort$start, Chr02_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(MvSl/Mint)",main="Chr02")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()

#autosome chromosome3
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMint_Chr03_ratio_20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr03_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr03_sort$ratio)),20)
plot(Chr03_sort$start, Chr03_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Log2(MvSl/Mint)",main="Chr03")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()

#autosome chromosome4
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater2n_Chr04_ratio_20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr04_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr04_sort$ratio)),20)
plot(Chr04_sort$start, Chr04_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Ratio(A1water/dikaryon)",main="Chr04")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()

#autosome chromosome5
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater2n_Chr05_ratio_20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr05_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr05_sort$ratio)),20)
plot(Chr05_sort$start, Chr05_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Ratio(A1water/dikaryon)",main="Chr05")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()

#autosome chromosome6
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidrich2n_Chr06_ratio_20genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr06_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr06_sort$ratio)),20)
plot(Chr06_sort$start, Chr06_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,4), xlab="Position(bp)", ylab="Ratio(A1water/dikaryon)",main="Chr06")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(h=1,lty=2,col="red")
abline(h=-1,lty=2,col="red")
abline(h=-2,lty=2,col="blue")
abline(h=2,lty=2,col="blue")
dev.off()



