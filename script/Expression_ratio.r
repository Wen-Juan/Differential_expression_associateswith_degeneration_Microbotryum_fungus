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

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amgonad_strict/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amgonad_strict/LogCPM_0.05_Amgonad copy.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

map.data <- subset(kdata, kdata$start!='NA')
chr1.1.data <- data.frame(subset(map.data, map.data$chr=='Chr01')) 
chr1.2.data <- data.frame(subset(map.data, map.data$chr=='Chr02')) 
chr1.data <- rbind(chr1.1.data,chr1.2.data)
chr2.1.data <- subset(map.data, map.data$chr!='Chr01') 
chr2.data<- subset(chr2.1.data, chr2.1.data$chr!='Chr02') 

#group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY male group at G46
#group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY male group at G46
#group1.expr.minus <- (map.data$Am1_434 + map.data$Am2_434 + map.data$Am4_434) / 3 #XY male group at G43
group1.expr.minus <- (map.data$A15MT1 + map.data$A17MT1 + map.data$A6MT1 + map.data$A8MT1 + map.data$A12MT1) / 5 #XY male group in gonad
#group1.expr.minus <- (map.data$A10MB + map.data$A15MB + map.data$A16MB + map.data$A17MB + map.data$A8MB) /5 #XY in brain
#group1.expr.minus <- (map.data$A12ML1 + map.data$A17ML1 + map.data$A15ML1 +map.data$A8ML1)/4 #XY linver

#group2.expr.minus <- (map.data$Am5_461) #Sex reversal XX male
#group2.expr.minus <- (map.data$Am2_464+map.data$Am6_464)/2 #XX female group at G46
#group2.expr.minus <- (map.data$Am2_433 + map.data$Am4_435 + map.data$Am5_433) / 3 #XX female group at G43
group2.expr.minus <- (map.data$A10FO1 + map.data$A17FO1 + map.data$A2FO2 + map.data$A8FO2 + map.data$A6FO1) / 5 #XY male group in gonad
#group2.expr.minus <- (map.data$A10FB + map.data$A12FB + map.data$A15FB + map.data$A16FB + map.data$A17FB) /5 #XX in brain
#group2.expr.minus <- (map.data$A12FL1 + map.data$A17FL1 + map.data$A2FL1 + map.data$A7FL1 +map.data$A8FL1)/5 #XX linver

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))


map.data$ratio <- log2(group1.expr/group2.expr)
#map.data$ratio <- group1.expr/group2.expr
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

#gene expression ratio Log2(XY/XX)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/log2ratio_XYXXsexreversal_10chr.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-6,6)) + 
  geom_boxplot(notch = TRUE) +
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


#gene expression ratio Log2(XY/XX)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/log2ratio_XYXXsexreversal_allauto.pdf", width=8, height=8)
ggplot(map.data, aes(x=chrone, y=ratio, fill=chrone)) +
  scale_fill_manual(values = c("grey","red","red")) +
  scale_y_continuous(limits = c(-0.5,0.5)) + 
  geom_boxplot(notch = TRUE,outlier.shape=NA,)+
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#stats
##
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 291820, p-value = 0.9756
withoutchr01 <- subset(map.data, map.data$chr!='Chr01')
withoutchr0102 <- subset(withoutchr01, withoutchr01$chr!='Chr02')

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],withoutchr01$ratio[withoutchr01$chr!='Chr02'],exact = FALSE) 
##W = 1525200, p-value = 0.513
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],withoutchr01$ratio[withoutchr01$chr!='Chr02'],exact = FALSE) 
##W = 1350200, p-value = 0.4937

###liver
#ggplot of boxplot
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/log2ratio_liver.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-8,8)) + 
  geom_boxplot() +
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###G46 no sex reversals

####
##
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 269820, p-value = 0.6143
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1734400, p-value = 0.9119
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1565500, p-value = 0.453
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 3512100, p-value = 0.4876

###
###G46 no sex reversals

####
##G43
##
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 293860, p-value = 0.2822
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1972500, p-value = 0.8631
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1796200, p-value = 0.1252
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 4026500, p-value = 0.2345He

#gonad

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 183190, p-value = 0.164
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1113200, p-value = 0.6068
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 951590, p-value = 0.1998
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 2140700, p-value = 0.5619

##brain

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 329790, p-value = 0.8941
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 2161000, p-value = 0.4588
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1877400, p-value = 0.5987
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 4266700, p-value = 0.3534

##liver

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 159060, p-value = 0.8283
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1036600, p-value = 0.7458
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 853180, p-value = 0.5964
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 2013500, p-value = 0.4988

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



dn_ds <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds/amm88perc_annotation_dnds.txt", header = T)
str(dn_ds)

##Assigning Transcripts to Chromosomes###
Chr1 <- subset(dn_ds, dn_ds$chr=="Chr01")
Chr2 <- subset(dn_ds, dn_ds$chr=="Chr02")
Chr3 <- subset(dn_ds, dn_ds$chr=="Chr03")
Chr4 <- subset(dn_ds, dn_ds$chr=="Chr04")
Chr5 <- subset(dn_ds, dn_ds$chr=="Chr05")
Chr6 <- subset(dn_ds, dn_ds$chr=="Chr06")
Chr7 <- subset(dn_ds, dn_ds$chr=="Chr07")
Chr8 <- subset(dn_ds, dn_ds$chr=="Chr08")
Chr9 <- subset(dn_ds, dn_ds$chr=="Chr09")
Chr10 <- subset(dn_ds, dn_ds$chr=="Chr10")

## Sort According to Position ####
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

###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

Chr1RT<- rollmean(smooth(Chr01_sort$dnds),40)
Chr2RT<- rollmean(smooth(Chr02_sort$dnds),40)
Chr3RT<- rollmean(smooth(Chr03_sort$dnds),40)
Chr4RT<- rollmean(smooth(Chr04_sort$dnds),40)
Chr5RT<- rollmean(smooth(Chr05_sort$dnds),40)
Chr6RT<- rollmean(smooth(Chr06_sort$dnds),40)
Chr7RT<- rollmean(smooth(Chr07_sort$dnds),40)
Chr8RT<- rollmean(smooth(Chr08_sort$dnds),40)
Chr9RT<- rollmean(smooth(Chr09_sort$dnds),40)
Chr10RT<- rollmean(smooth(Chr10_sort$dnds),40)


Rt<- c(Chr3RT,Chr4RT,Chr5RT,Chr6RT,Chr7RT,Chr8RT,Chr9RT,Chr10RT)


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

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/chr1_dnds.pdf", width=7,height=5)

Chr_pos <- rollmean(smooth(Chr01_sort$start),40)
Chr_dnds <- rollmean(smooth(Chr01_sort$dnds),40)
plot(Chr01_sort$start, Chr01_sort$dnds,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(0,0.4), xlab="Position(bp)", ylab="dN/dS",main="Chr01")
lines(Chr_pos, Chr_dnds,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
abline(v=116443467, col="blue",lwd=2,lty=3)
dev.off()

###gene expression at G46
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

## Sort According to Position ####
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

###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

Chr1RT<- rollapply(Chr01_sort$ratio, 40, mean, na.rm = TRUE) #need to switch to rollapply, cause rollmean does not handle NA values.
Chr2RT<- rollapply(Chr02_sort$ratio, 40, mean, na.rm = TRUE)
Chr3RT<- rollapply(Chr03_sort$ratio, 40, mean, na.rm = TRUE)
Chr4RT<- rollapply(Chr04_sort$ratio, 40, mean, na.rm = TRUE)
Chr5RT<- rollapply(Chr05_sort$ratio, 40, mean, na.rm = TRUE)
Chr6RT<- rollapply(Chr06_sort$ratio, 40, mean, na.rm = TRUE)
Chr7RT<- rollapply(Chr07_sort$ratio, 40, mean, na.rm = TRUE)
Chr8RT<- rollapply(Chr08_sort$ratio, 40, mean, na.rm = TRUE)
Chr9RT<- rollapply(Chr09_sort$ratio, 40, mean, na.rm = TRUE)
Chr10RT<- rollapply(Chr10_sort$ratio, 40, mean, na.rm = TRUE)


Rt<- c(Chr3RT,Chr4RT,Chr5RT,Chr6RT,Chr7RT,Chr8RT,Chr9RT,Chr10RT)


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

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/chr10_expratio_G46.pdf", width=7,height=5)

Chr_pos <- rollapply(Chr10_sort$start, 40, mean, na.rm = TRUE)
Chr_ratio <- rollapply(Chr10_sort$ratio, 40, mean, na.rm = TRUE)
plot(Chr10_sort$start, Chr10_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-6,6), xlab="Position(bp)", ylab="Log2(male:female)",main="Chr10")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
#abline(v=116443467, col="blue",lwd=2,lty=3)
dev.off()

#gene expression ratio for gonad
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

## Sort According to Position ####
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

Chr1RT<- rollapply(Chr01_sort$ratio, 40, mean, na.rm = TRUE) #need to switch to rollapply, cause rollmean does not handle NA values.
Chr2RT<- rollapply(Chr02_sort$ratio, 40, mean, na.rm = TRUE)
Chr3RT<- rollapply(Chr03_sort$ratio, 40, mean, na.rm = TRUE)
Chr4RT<- rollapply(Chr04_sort$ratio, 40, mean, na.rm = TRUE)
Chr5RT<- rollapply(Chr05_sort$ratio, 40, mean, na.rm = TRUE)
Chr6RT<- rollapply(Chr06_sort$ratio, 40, mean, na.rm = TRUE)
Chr7RT<- rollapply(Chr07_sort$ratio, 40, mean, na.rm = TRUE)
Chr8RT<- rollapply(Chr08_sort$ratio, 40, mean, na.rm = TRUE)
Chr9RT<- rollapply(Chr09_sort$ratio, 40, mean, na.rm = TRUE)
Chr10RT<- rollapply(Chr10_sort$ratio, 40, mean, na.rm = TRUE)


Rt<- c(Chr3RT,Chr4RT,Chr5RT,Chr6RT,Chr7RT,Chr8RT,Chr9RT,Chr10RT)


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

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/chr10_expratio_gonad.pdf", width=7,height=5)

Chr_pos <- rollapply(Chr10_sort$start, 40, mean, na.rm = TRUE)
Chr_ratio <- rollapply(Chr10_sort$ratio, 40, mean, na.rm = TRUE)
plot(Chr10_sort$start, Chr10_sort$ratio,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-10,10), xlab="Position(bp)", ylab="Log2(male:female)",main="Chr10")
lines(Chr_pos, Chr_ratio,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
#abline(v=116443467, col="blue",lwd=2,lty=3)

dev.off()


##WJ original code without confidential internal
map.data <- dn_ds
chr1_1_start <-map.data$start[map.data$chr=='Chr01']
chr_1_end <- max(map.data$start[map.data$chr=='Chr01'])

zoo.dat <- zoo(map.data$ratio[map.data$chr=='Chr01'], c(chr1_1_start,chr_1_end))
y <- rollapply(zoo.dat, 40, FUN = mean, align = 'center', na.rm=TRUE) 


pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_chr01_malesexreveral_ratio.pdf", width=8, height=8)
par(mar=c(5,5,4,3)+0.6) 
#plot(c(0, chr_1_end), c(0,2), axes=F, lwd=2, xlab="Position (bp)", ylab="gene expression ratio log2(XY0/XX)", cex.axis=1.5, cex.lab=1.2, col="white")
plot(c(0, chr_1_end), c(-6,6), axes=F, lwd=2, xlab="Position (bp)", ylab="gene expression ratio Log2(XY/XX)", cex.axis=1.5, cex.lab=1.2, col="white")
axis(1,c(0, 30000000,60000000,90000000, 120000000, 150000000, 180000000, 200000000))
axis(2, c(-6,-4,-2,0,2,4,6))
points(map.data$start[map.data$chr=='Chr01'], map.data$ratio[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="gray50", main="",cex.axis=1.5)
#points(DC.data1dmrt$start[DC.data1dmrt$chr=='Chr01'],DC.data1dmrt$dN_dS[DC.data1dmrt$chr=='Chr01'], pch=20, lwd=3, type="p",col="red", main="",cex.axis=1.5)
lines(y,col="blue",lwd=4)
abline(v=116443467, col="blue",lwd=2,lty=3)
dev.off()

