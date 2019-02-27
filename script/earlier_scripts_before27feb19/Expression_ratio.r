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

datapath <- '/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/haploidwater_v1/A1allgenes/'
kdata <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/A1allgenes/LogCPM_MintMvsl_singlecopy_genegain_all.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- data.frame(subset(map.data, map.data$chr=='aMAT')) 
chr2.data <- subset(map.data, map.data$chr!='aMAT') 

group1.expr.minus <- (map.data$water1 + map.data$water2) / 2
group2.expr.minus <- (map.data$di1 + map.data$di2) / 2

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- log2(group1.expr/group2.expr)
#map.data$ratio <- group1.expr/group2.expr
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mint_Mvsl_singleortholog_a1a2_expression_ratio.pdf", width=8, height=8)


write.table(map.data, file = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mint_mvsl_singlecopyortholog_expratio.txt", sep='\t', row.names = TRUE, col.names = TRUE)
#gene expression ratio Log2(XY/XX)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_log2ratio_A1_A2.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red")) +
  scale_y_continuous(limits = c(-0.3,0.3)) + 
  geom_boxplot(notch = FALSE,outlier.shape=NA,)+
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of A1/A2 gene expression') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/expression_MintMvsl_singlecopy_genegain_haploidvsdikaryon.pdf", width=8, height=8)
ggplot(map.data, aes(x = chr, y=logFC.water.diploid, fill=ortholog)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), name="Ortholog or not",labels=c("gene gain","single ortholog")) +
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0.7), width=0.7,alpha=0.8) +
  coord_cartesian(ylim = c(-3,-2,-1,0,1,2,3)) +
  scale_x_discrete(labels=c("MAT","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","rest")) +
  geom_vline(aes(xintercept=10), linetype="blank") +
  labs(x='Chromosome', y='Log2(haploid/dikaryon)')
dev.off()

#stats
##
wilcox.test(map.data$ratio[map.data$chr=='aMAT'],map.data$ratio[map.data$chr!='aMAT'],exact = FALSE) 
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

###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

MATRT<- rollmean(smooth(na.approx(MAT_sort$ratio)),10)
Chr1RT<- rollmean(smooth(na.approx(Chr01_sort$ratio)),10)
Chr2RT<- rollmean(smooth(na.approx(Chr02_sort$ratio)),10)
Chr3RT<- rollmean(smooth(na.approx(Chr03_sort$ratio)),10)
Chr4RT<- rollmean(smooth(na.approx(Chr04_sort$ratio)),10)
Chr5RT<- rollmean(smooth(na.approx(Chr05_sort$ratio)),10)
Chr6RT<- rollmean(smooth(na.approx(Chr06_sort$ratio)),10)
Chr7RT<- rollmean(smooth(na.approx(Chr07_sort$ratio)),10)
Chr8RT<- rollmean(smooth(na.approx(Chr08_sort$ratio)),10)
Chr9RT<- rollmean(smooth(na.approx(Chr09_sort$ratio)),10)
Chr10RT<- rollmean(smooth(na.approx(Chr10_sort$ratio)),10)
Chr11RT<- rollmean(smooth(na.approx(Chr11_sort$ratio)),10)
Chr12RT<- rollmean(smooth(na.approx(Chr12_sort$ratio)),10)
Chr13RT<- rollmean(smooth(na.approx(Chr13_sort$ratio)),10)
Chr14RT<- rollmean(smooth(na.approx(Chr14_sort$ratio)),10)
Chr15RT<- rollmean(smooth(na.approx(Chr15_sort$ratio)),10)
Chr16RT<- rollmean(smooth(na.approx(Chr16_sort$ratio)),10)
Chr17RT<- rollmean(smooth(na.approx(Chr17_sort$ratio)),10)

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

RMpalette <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")

#mating type chromosome
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mint_Mvsl_genegainsMAT_expratio.pdf", width=7,height=5)

MAT_ortholog <- subset(MAT_sort, MAT_sort$ortholog != "genegain")
Chr_pos <- rollmean(smooth(MAT_ortholog$start),10)
Chr_ratio <- rollmean(smooth(na.approx(MAT_ortholog$ratio)),10)
plot(MAT_ortholog$start, MAT_ortholog$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,5), xlab="Position(bp)", ylab="Log2(haploidwater/dikaryon)",main="Mint and MvSl gene gains on MAT")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
#abline(h=1,lty=2,col="red")
#abline(h=-1,lty=2,col="red")

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
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mint_Mvsl_genegainsChr01_expratio.pdf", width=7,height=5)

chr01_ortholog <- subset(Chr01_sort, Chr01_sort$ortholog == "genegain")
Chr_pos <- rollmean(smooth(chr01_ortholog$start),10)
Chr_ratio <- rollmean(smooth(na.approx(chr01_ortholog$ratio)),10)
plot(chr01_ortholog$start, chr01_ortholog$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-4,5), xlab="Position(bp)", ylab="Log2(haploidwater/dikaryon)",main="Mint and MvSl gene gains on Chr01")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

##autosome chromosome2
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_Chr02_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr02_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr02_sort$ratio)),20)
plot(Chr02_sort$start, Chr02_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1.5,1.5), xlab="Position(bp)", ylab="Log2-ratio(A1/A2)",main="Chr02")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

#autosome chromosome3
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_Chr03_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr03_sort$start),40)
Chr_ratio <- rollmean(smooth(na.approx(Chr03_sort$ratio)),40)
plot(Chr03_sort$start, Chr03_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-2.5,2.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="Chr03")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

#autosome chromosome4
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_Chr04_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr04_sort$start),40)
Chr_ratio <- rollmean(smooth(na.approx(Chr04_sort$ratio)),40)
plot(Chr04_sort$start, Chr04_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1.5,1.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="Chr04")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

#autosome chromosome5
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_Chr05_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr05_sort$start),40)
Chr_ratio <- rollmean(smooth(na.approx(Chr05_sort$ratio)),40)
plot(Chr05_sort$start, Chr05_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1.5,1.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="Chr05")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

#autosome chromosome6
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/haploidwater_Chr06_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr06_sort$start),40)
Chr_ratio <- rollmean(smooth(na.approx(Chr06_sort$ratio)),40)
plot(Chr06_sort$start, Chr06_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-2.5,2.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="Chr06")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()


