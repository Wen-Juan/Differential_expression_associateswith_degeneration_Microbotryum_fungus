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

datapath <- '/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/haploidwater_new5/'
kdata <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/haploidwater_new5/LogCPM_0.05_haploidwater.txt",header = T)
str(kdata)

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- data.frame(subset(map.data, map.data$chr=='aMAT')) 
chr2.data <- subset(map.data, map.data$chr!='aMAT') 

group1.expr.minus <- (map.data$A1_1 + map.data$A1_2) / 2
group2.expr.minus <- (map.data$A2_1 + map.data$A2_2) / 2

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- abs(log2(group1.expr/group2.expr))
#map.data$ratio <- group1.expr/group2.expr
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

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

##Assigning Transcripts to Chromosomes###
MAT <- subset(map.data, map.data$chr=="aMAT")
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
Chr18 <- subset(map.data, map.data$chr=="Chr18")


## Sort According to Position ####
MAT_sort <- MAT[order(MAT$start),]
Chr01_sort <- Chr1[order(Chr1$start),] 
Chr02_sort <- Chr2[order(Chr2$start),] 
Chr03_sort <- Chr3[order(Chr3$start),] 
Chr04_sort <- Chr4[order(Chr4$start),] 
Chr05_sort <- Chr5[order(Chr5$start),] 
Chr06_sort <- Chr6[order(Chr6$start),] 
Chr07_sort <- Chr7[order(Chr7$start),] 
Chr08_sort <- Chr8[order(Chr8$start),] 
Chr09_sort <- Chr9[order(Chr9$start),] 
Chr10_sort <- Chr10[order(Chr10$start),] 
Chr11_sort <- Chr11[order(Chr11$start),] 
Chr12_sort <- Chr12[order(Chr12$start),] 
Chr13_sort <- Chr13[order(Chr13$start),] 
Chr14_sort <- Chr6[order(Chr14$start),] 
Chr15_sort <- Chr6[order(Chr15$start),] 
Chr16_sort <- Chr6[order(Chr16$start),] 
Chr17_sort <- Chr6[order(Chr17$start),] 

###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

MATRT<- rollmean(smooth(na.approx(MAT_sort$ratio)),20)
Chr1RT<- rollmean(smooth(na.approx(Chr01_sort$ratio)),20)
Chr2RT<- rollmean(smooth(na.approx(Chr02_sort$ratio)),20)
Chr3RT<- rollmean(smooth(na.approx(Chr03_sort$ratio)),20)
Chr4RT<- rollmean(smooth(na.approx(Chr04_sort$ratio)),20)
Chr5RT<- rollmean(smooth(na.approx(Chr05_sort$ratio)),20)
Chr6RT<- rollmean(smooth(na.approx(Chr06_sort$ratio)),20)
Chr7RT<- rollmean(smooth(na.approx(Chr07_sort$ratio)),20)
Chr8RT<- rollmean(smooth(na.approx(Chr08_sort$ratio)),20)
Chr9RT<- rollmean(smooth(na.approx(Chr09_sort$ratio)),20)
Chr10RT<- rollmean(smooth(na.approx(Chr10_sort$ratio)),20)
Chr11RT<- rollmean(smooth(na.approx(Chr11_sort$ratio)),20)
Chr12RT<- rollmean(smooth(na.approx(Chr12_sort$ratio)),20)
Chr13RT<- rollmean(smooth(na.approx(Chr13_sort$ratio)),20)
Chr14RT<- rollmean(smooth(na.approx(Chr14_sort$ratio)),20)
Chr15RT<- rollmean(smooth(na.approx(Chr15_sort$ratio)),20)
Chr16RT<- rollmean(smooth(na.approx(Chr16_sort$ratio)),20)
Chr17RT<- rollmean(smooth(na.approx(Chr17_sort$ratio)),20)

Rt<- c(MATRT,Chr1RT,Chr2RT,Chr3RT,Chr4RT,Chr5RT,Chr6RT,Chr7RT,Chr8RT,Chr9RT,Chr10RT,Chr11RT,Chr12RT,Chr13RT,Chr14RT,Chr15RT,Chr16RT,Chr17RT)

myfunction <- function(i){
  Info <- sample(i,1,replace=FALSE)
  return(Info)
}

my.perm <- c()
for(i in 1:10^3){ my.perm[i] <- myfunction(Rt) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]

RMpalette <- c("#f0f9e8", "#bae4bc", "black", "#43a2ca", "red")

#mating type chromosome
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MAT_nn_expressionratio.pdf", width=7,height=5)

MAT <- subset(MAT_sort, MAT_sort$chr == "aMAT")
Chr_pos <- rollmean(smooth(MAT$start),20)
Chr_ratio <- rollmean(smooth(na.approx(MAT$logFC.A1.A2)),20)
plot(MAT$start, MAT$logFC.A1.A2,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-10,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on MAT")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
#abline(h=lowCI,lty=2)
#abline(h=highCI,lty=2)
#abline(h=1,lty=2,col="red")
#abline(h=-1,lty=2,col="red")

dev.off()

#MAT absolute log2
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Abs_MAT_nn_expressionratio.pdf", width=7,height=5)

MAT <- subset(MAT_sort, MAT_sort$chr == "aMAT")
Chr_pos <- rollmean(smooth(MAT$start),20)
Chr_ratio <- rollmean(smooth(na.approx(MAT$absLogFC )),20)
plot(MAT$start, MAT$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, xlim=c(0,3500000), ylim=c(0,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on MAT")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

wilcox.test(MAT$ratio[MAT$ortholog=='genegain'],MAT$ratio[MAT$ortholog=='sortholog'],exact = FALSE) 
#W = 127120, p-value < 2.2e-16
wilcox.test(Chr1$ratio[Chr1$ortholog=='genegain'],Chr1$ratio[Chr1$ortholog=='sortholog'],exact = FALSE) 
#W = 123760, p-value = 0.004411
wilcox.test(Chr2$ratio[Chr2$ortholog=='genegain'],Chr2$ratio[Chr2$ortholog=='sortholog'],exact = FALSE) 
#W = 62361, p-value = 1.176e-07
wilcox.test(Chr3$ratio[Chr3$ortholog=='genegain'],Chr3$ratio[Chr3$ortholog=='sortholog'],exact = FALSE) 
#W = 54221, p-value = 0.0001812
wilcox.test(Chr4$ratio[Chr4$ortholog=='genegain'],Chr4$ratio[Chr4$ortholog=='sortholog'],exact = FALSE) 
#W = 36785, p-value = 0.01385
wilcox.test(Chr5$ratio[Chr5$ortholog=='genegain'],Chr5$ratio[Chr5$ortholog=='sortholog'],exact = FALSE) 
#W = 29252, p-value = 0.09072
wilcox.test(Chr6$ratio[Chr6$ortholog=='genegain'],Chr6$ratio[Chr6$ortholog=='sortholog'],exact = FALSE) 
#W = 26027, p-value = 0.005781
wilcox.test(Chr7$ratio[Chr7$ortholog=='genegain'],Chr7$ratio[Chr7$ortholog=='sortholog'],exact = FALSE) 
#W = 26027, p-value = 0.005781
wilcox.test(Chr8$ratio[Chr8$ortholog=='genegain'],Chr8$ratio[Chr8$ortholog=='sortholog'],exact = FALSE) 
#W = 26733, p-value = 0.6509
wilcox.test(Chr9$ratio[Chr9$ortholog=='genegain'],Chr9$ratio[Chr9$ortholog=='sortholog'],exact = FALSE) 
#W = 22003, p-value = 0.0002844
wilcox.test(Chr10$ratio[Chr10$ortholog=='genegain'],Chr10$ratio[Chr10$ortholog=='sortholog'],exact = FALSE) 
#W = 19022, p-value = 0.08431
wilcox.test(Chr11$ratio[Chr11$ortholog=='genegain'],Chr11$ratio[Chr11$ortholog=='sortholog'],exact = FALSE) 
#W = 18539, p-value = 0.02587
wilcox.test(Chr12$ratio[Chr12$ortholog=='genegain'],Chr12$ratio[Chr12$ortholog=='sortholog'],exact = FALSE)
#W = 20210, p-value = 0.7273
wilcox.test(Chr13$ratio[Chr13$ortholog=='genegain'],Chr13$ratio[Chr13$ortholog=='sortholog'],exact = FALSE)
#W = 11970, p-value = 0.4994
wilcox.test(Chr14$ratio[Chr14$ortholog=='genegain'],Chr14$ratio[Chr14$ortholog=='sortholog'],exact = FALSE)
#W = 3768, p-value = 0.08985
wilcox.test(Chr15$ratio[Chr15$ortholog=='genegain'],Chr15$ratio[Chr15$ortholog=='sortholog'],exact = FALSE)
#W = 675, p-value = 0.001431
wilcox.test(Chr16$ratio[Chr16$ortholog=='genegain'],Chr16$ratio[Chr16$ortholog=='sortholog'],exact = FALSE)
#W = 3, p-value = 0.04296
wilcox.test(Chr17$ratio[Chr17$ortholog=='genegain'],Chr17$ratio[Chr17$ortholog=='sortholog'],exact = FALSE)
#W = 17, p-value = 0.1706
wilcox.test(Chr18$ratio[Chr18$ortholog=='genegain'],Chr18$ratio[Chr18$ortholog=='sortholog'],exact = FALSE)
#W = 17, p-value = 0.1706

#autosome chromosome1
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Chr01_nn_expressionratio.pdf", width=7,height=5)

Chr01 <- subset(Chr01_sort, Chr01_sort$chr == "Chr01")
Chr_pos <- rollmean(smooth(Chr01$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr01$logFC.A1.A2)),20)
plot(Chr01$start, Chr01$logFC.A1.A2,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-10,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr01")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])

dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/asb_Chr01_nn_expressionratio.pdf", width=7,height=5)

Chr01 <- subset(Chr01_sort, Chr01_sort$chr == "Chr01")
Chr_pos <- rollmean(smooth(Chr01$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr01$absLogFC)),20)
plot(Chr01$start, Chr01$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, xlim=c(0,3500000),ylim=c(0,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr01")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

##autosome chromosome2
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr02_nn_expressionratio.pdf", width=7,height=5)

Chr02 <- subset(Chr02_sort, Chr02_sort$chr == "Chr02")
Chr_pos <- rollmean(smooth(Chr02$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr02$absLogFC)),20)
plot(Chr02$start, Chr02$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10), xlab="Position(bp)",  xlim=c(0,3500000),ylab="Log2(a1/a2)",main="Expression ratio on Chr02")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#autosome chromosome3
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr03_nn_expressionratio.pdf", width=7,height=5)

Chr03 <- subset(Chr03_sort, Chr03_sort$chr == "Chr03")
Chr_pos <- rollmean(smooth(Chr03$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr03$absLogFC)),20)
plot(Chr03$start, Chr03$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr03")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#autosome chromosome4
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Chr04_nn_expressionratio.pdf", width=7,height=5)

Chr04 <- subset(Chr04_sort, Chr04_sort$chr == "Chr04")
Chr_pos <- rollmean(smooth(Chr04$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr04$logFC.A1.A2)),20)
plot(Chr04$start, Chr04$logFC.A1.A2,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-10,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr04")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/asb_Chr04_nn_expressionratio.pdf", width=7,height=5)

Chr04 <- subset(Chr04_sort, Chr04_sort$chr == "Chr04")
Chr_pos <- rollmean(smooth(Chr04$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr04$absLogFC)),20)
plot(Chr04$start, Chr04$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr04")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()
#autosome chromosome5
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr05_nn_expressionratio.pdf", width=7,height=5)

Chr05 <- subset(Chr05_sort, Chr05_sort$chr == "Chr05")
Chr_pos <- rollmean(smooth(Chr05$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr05$absLogFC)),20)
plot(Chr05$start, Chr05$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr05")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#autosome chromosome6
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr06_nn_expressionratio.pdf", width=7,height=5)

Chr06 <- subset(Chr06_sort, Chr06_sort$chr == "Chr06")
Chr_pos <- rollmean(smooth(Chr06$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr06$absLogFC)),20)
plot(Chr06$start, Chr06$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr06")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#auto chr07
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr07_nn_expressionratio.pdf", width=7,height=5)

Chr07 <- subset(Chr07_sort, Chr07_sort$chr == "Chr07")
Chr_pos <- rollmean(smooth(Chr07$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr07$absLogFC)),20)
plot(Chr07$start, Chr07$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr07")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#auto chr08
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr08_nn_expressionratio.pdf", width=7,height=5)

Chr08 <- subset(Chr08_sort, Chr08_sort$chr == "Chr08")
Chr_pos <- rollmean(smooth(Chr08$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr08$absLogFC)),20)
plot(Chr08$start, Chr08$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr08")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#chr09
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr09_nn_expressionratio.pdf", width=7,height=5)

Chr09 <- subset(Chr09_sort, Chr09_sort$chr == "Chr09")
Chr_pos <- rollmean(smooth(Chr09$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr09$absLogFC)),20)
plot(Chr09$start, Chr09$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr09")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#chr10
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr10_nn_expressionratio.pdf", width=7,height=5)

Chr10 <- subset(Chr10_sort, Chr10_sort$chr == "Chr10")
Chr_pos <- rollmean(smooth(Chr10$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr10$absLogFC)),20)
plot(Chr10$start, Chr10$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr10")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#chr11
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr11_nn_expressionratio.pdf", width=7,height=5)

Chr11 <- subset(Chr11_sort, Chr11_sort$chr == "Chr11")
Chr_pos <- rollmean(smooth(Chr11$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr11$absLogFC)),20)
plot(Chr11$start, Chr11$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,10),  xlim=c(0,3500000),xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr11")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#chr12
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr12_nn_expressionratio.pdf", width=7,height=5)

Chr12 <- subset(Chr12_sort, Chr12_sort$chr == "Chr12")
Chr_pos <- rollmean(smooth(Chr12$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr12$absLogFC)),20)
plot(Chr12$start, Chr12$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20, xlim=c(0,3500000), ylim=c(0,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr12")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()

#chr13
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/abs_Chr13_nn_expressionratio.pdf", width=7,height=5)

Chr13 <- subset(Chr13_sort, Chr13_sort$chr == "Chr13")
Chr_pos <- rollmean(smooth(Chr13$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr13$absLogFC)),20)
plot(Chr13$start, Chr13$absLogFC,col=alpha(RMpalette[3], 0.5),pch=20,  xlim=c(0,3500000), ylim=c(0,10), xlab="Position(bp)", ylab="Log2(a1/a2)",main="Expression ratio on Chr13")
lines(Chr_pos, Chr_ratio,type="l",lwd=4, col=RMpalette[5])

dev.off()
