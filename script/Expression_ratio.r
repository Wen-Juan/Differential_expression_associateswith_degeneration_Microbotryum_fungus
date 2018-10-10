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

datapath <- '/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/haploidrich/'
kdata <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/haploidrich/LogCPM_0.05_haploidrich copy.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- data.frame(subset(map.data, map.data$chr=='aMAT')) 
chr2.data <- subset(map.data, map.data$chr!='aMAT') 

group1.expr.minus <- (map.data$A1_1 + map.data$A1_2) / 2
group2.expr.minus <- (map.data$A2_1 + map.data$A2_2) / 2

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- log2(group1.expr/group2.expr)
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

#gene expression ratio Log2(XY/XX)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/log2ratio_A1_A2.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-0.3,0.3)) + 
  geom_boxplot(notch = TRUE,outlier.shape=NA,)+
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of A1/A2 gene expression') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#stats
##
wilcox.test(map.data$ratio[map.data$chr=='aMAT'],map.data$ratio[map.data$chr!='aMAT'],exact = FALSE) 
##W = 1260900, p-value = 0.9041

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
aMAT <- subset(map.data, map.data$chr=="aMAT")
Chr1 <- subset(map.data, map.data$chr=="Chr01")
Chr2 <- subset(map.data, map.data$chr=="Chr02")
Chr3 <- subset(map.data, map.data$chr=="Chr03")
Chr4 <- subset(map.data, map.data$chr=="Chr04")
Chr5 <- subset(map.data, map.data$chr=="Chr05")
Chr6 <- subset(map.data, map.data$chr=="Chr06")


## Sort According to Position ####
aMAT_sort <- aMAT[order(aMAT$start),]
Chr01_sort <- Chr1[order(Chr1$start),] 
Chr02_sort <- Chr2[order(Chr2$start),] 
Chr03_sort <- Chr3[order(Chr3$start),] 
Chr04_sort <- Chr4[order(Chr4$start),] 
Chr05_sort <- Chr5[order(Chr5$start),] 
Chr06_sort <- Chr6[order(Chr6$start),] 


###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

aMATRT<- rollmean(smooth(na.approx(aMAT_sort$ratio)),40)
Chr1RT<- rollmean(smooth(na.approx(Chr01_sort$ratio)),40)
Chr2RT<- rollmean(smooth(na.approx(Chr02_sort$ratio)),40)
Chr3RT<- rollmean(smooth(na.approx(Chr03_sort$ratio)),40)
Chr4RT<- rollmean(smooth(na.approx(Chr04_sort$ratio)),40)
Chr5RT<- rollmean(smooth(na.approx(Chr05_sort$ratio)),40)
Chr6RT<- rollmean(smooth(na.approx(Chr06_sort$ratio)),40)

Rt<- c(aMATRT,Chr1RT,Chr2RT,Chr3RT,Chr4RT,Chr5RT,Chr6RT)

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
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MAT_ratio.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(aMAT_sort$start),10)
Chr_ratio <- rollmean(smooth(na.approx(aMAT_sort$ratio)),10)
plot(aMAT_sort$start, aMAT_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1.5,1.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="MAT")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

#autosome chromosome1
pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Chr01_ratio_sliding40genes.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr01_sort$start),20)
Chr_ratio <- rollmean(smooth(na.approx(Chr01_sort$ratio)),20)
plot(Chr01_sort$start, Chr01_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-1.5,1.5), xlab="Position(bp)", ylab="Ratio(A1/A2)",main="Chr01")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

