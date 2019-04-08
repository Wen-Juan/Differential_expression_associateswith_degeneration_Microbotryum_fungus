# install packages and load R libraries 
install.packages('ggfortify')
library('ggfortify')

install.packages("missMDA")
library(missMDA)

source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

install.packages("lmerTest")
library(lme4)
library(lmerTest)

install.packages("statmod")
library('statmod')

install.packages("lattice")
library(lattice)
install.packages("Amelia")
library('Amelia')
install.packages("mice")
library(mice)


#load dataset
PCA_5degen_all <- read.table ("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/PCA/Mvsl_a1a2_exp_compartment_all5traits_nonoriented.txt", header = TRUE)
str(PCA_5degen_all)
PCA_5degen_all <- mice(PCA_5degen_all, m=5, method='norm.boot')

##stats with log2 transformation
qqnorm(log10(PCA_5degen_all$abs), pch = 1, frame = FALSE)
qqline(log10(PCA_5degen_all$abs), col = "steelblue")
hist(log2(PCA_5degen_all$abs))

hist(PCA_5degen_all$abs)
cor.test(log2(PCA_5degen_all$abs), log2(PCA_5degen_all$totalTE+0.00000000001), met0hod = "spearm", alternative = "g")

plot(log2(PCA_5degen_all$abs) ~ PCA_5degen_all$totalTE)

###stats without data transformation for the responsible variable.
hist(sqrt(sqrt(PCA_5degen_all$abs)))
qqnorm(sqrt(sqrt(PCA_5degen_all$abs)), pch = 1, frame = FALSE)

y1 <- glm(sqrt(sqrt(sqrt(abs))) ~ totalTE*youngold + diffGC3*youngold +dn*youngold +ratioprot*youngold +intronnrdiff*youngold+ intronmeandiff*youngold, data =PCA_5degen_all, family = gaussian)
summary(y1)

y2 <- glm(abs ~ totalTE*youngold + diffGC3*youngold +dn*youngold +ratioprot*youngold +intronnrdiff*youngold+ intronmeandiff*youngold, data =PCA_5degen_all, family=tweedie(var.power=1.8,link.power=0))
summary(y2)

######
######

##using variables of degeneration traits as column names
PCA_5degen1 <- as.data.frame(scale(PCA_5degen_all[,7:25]))
str(PCA_5degen1)
sum(is.na(PCA_5degen1))

pc <- pca(PCA_5degen1, nPcs=5, method="ppca")
imputed <- completeObs(pc)

####
pca <- prcomp(imputed, center = TRUE, scale. =  TRUE)
summary(pca)

pca$gencomp <- PCA_5degen_all[,2]
pca$de2 <- PCA_5degen_all[,6]
str(pca)

ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,  labels=rownames(mtcars), groups=mtcars.country) +
  scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue"))+
  ggtitle("PCA of mtcars dataset")+
  theme_minimal()+
  theme(legend.position = "bottom")

#pca1 and pca2
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca1pca2.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

#pca3 and pca4
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca3pca4.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, choices=c(3,4), obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

#pca5 and pca6
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_pca_pca5pca6.pdf", width=8, height=8)
ggbiplot(pca, ellipse = TRUE, circle=TRUE, choices=c(5,6), obs.scale = 1, var.scale = 1, labels = colnames(pca), groups = pca$gencomp) +
  scale_colour_manual(name="Genomic compartment", values= c("dark grey", "orange", "brown","green", "blue"))+
  theme_minimal()+
  theme(legend.position = "bottom")
dev.off()

