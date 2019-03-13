# install packages and load R libraries 
install.packages('ggfortify')
library('ggfortify')

#load dataset
PCA_5degen <- read.table ("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/PCA/Mvsl_common_all5degenerationtypes_sort.txt", header = TRUE)
str(PCA_5degen)

##using variables of degeneration traits as rownames
PCA_5degen1 <- data.frame(PCA_5degen[,7:25])
str(PCA_5degen1)
PCA_5degen2 <- data.frame(t(PCA_5degen1))
str(PCA_5degen2)
pcaData <- PCA_5degen2
pca <- prcomp(pcaData, scale. = TRUE)
summary(pca)

rownames(PCA_5degen2) <- c("genediff", "upk5diff", "upk10diff", "upk15diff", "upk20diff", "down5kdiff", "down10kdiff", "down15kdiff", "down20kdiff", "ratioprot", "ratiocds", "intron_nr_diff", "intron_mean_diff", "intron_total_diff", "GC0diff", "GC3diff","dn", "ds", "dnds")
str(PCA_5degen2)
PCA_5degen2$types <- c("TE", "TEup", "TEup", "TEup", "TEup", "TEdown", "TEdown", "TEdown", "TEdown", "Prot", "Prot", "intron", "intron", "intron", "GC0", "GC3","dn", "ds", "dnds") 
PCA_5degen2$type5 <- c("TE", "TE", "TE", "TE", "TE", "TE", "TE", "TE", "TE", "Prot", "Prot", "intron", "intron", "intron", "GC", "GC","dnds", "dnds", "dnds") 

PCA_5degen2_i <- data.frame(pca$x, types = PCA_5degen2$types, type5 = PCA_5degen2$type5)
str(PCA_5degen2_i)

autoplot(pca, data = PCA_5degen2_i, x=1, y=2, color = 'type5', shape = 'type5', label.size = 2,size=2, ylim=c(-1,1), xlim=c(-1,1))
autoplot(pca, data = PCA_5degen2_i, x=3, y=4,color = 'type5', shape = 'type5', label.size = 2,size=4, ylim=c(-2,2), xlim=c(-2,2))
autoplot(pca, data = PCA_5degen2_i, x=5, y=6,color = 'type5', shape = 'type5', label.size = 2,size=4, ylim=c(-2,2), xlim=c(-2,2))

##using variables of degeneration traits as column names
PCA_5degen1 <- data.frame(PCA_5degen[,7:25])
str(PCA_5degen1)
pcaData <- PCA_5degen1
pca <- prcomp(pcaData, scale. = TRUE)
summary(pca)

rownames(PCA_5degen2) <- c("genediff", "upk5diff", "upk10diff", "upk15diff", "upk20diff", "down5kdiff", "down10kdiff", "down15kdiff", "down20kdiff", "ratioprot", "ratiocds", "intron_nr_diff", "intron_mean_diff", "intron_total_diff", "GC0diff", "GC3diff","dn", "ds", "dnds")
str(PCA_5degen2)
PCA_5degen2$types <- c("TE", "TEup", "TEup", "TEup", "TEup", "TEdown", "TEdown", "TEdown", "TEdown", "Prot", "Prot", "intron", "intron", "intron", "GC0", "GC3","dn", "ds", "dnds") 
PCA_5degen2$type5 <- c("TE", "TE", "TE", "TE", "TE", "TE", "TE", "TE", "TE", "Prot", "Prot", "intron", "intron", "intron", "GC", "GC","dnds", "dnds", "dnds") 


#dev.off()
#pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/alltissues_pc3pc4.pdf")
#autoplot(pca, data = count_i, x=3, y=4, colour = 'stage', shape='sex', size=4, ylim=c(-0.3,0.3), xlim=c(-0.3,0.5))
#dev.off()

#pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/alltissues_pc5pc6.pdf")
#autoplot(pca, data = count_i, x=5, y=6, colour = 'stage', shape='sex', size=4, ylim=c(-0.3,0.5), xlim=c(-0.3,0.5))
#dev.off()

#ggplot (count_i, aes(x=PC1,y=PC2,col=stage,shape=sex)) + #this does not give variance percentage.
#  geom_point(size=3,alpha=0.8)+
#  scale_color_manual(values = c("orange", "red","pink","darkblue","blue", "brown", "green",  "grey"))+ 
#  theme_set(theme_bw(base_size=12)) +
#  theme(legend.justification=c(1,0), legend.position=c(0.95,0.05)) 


#ggplot (count_i, aes(x=PC3,y=PC4,col=stage,shape=sex,variance_percentage = TRUE)) + #this does not give variance percentage.
# geom_point(size=3,alpha=0.8) +
#  scale_color_manual(values = c("orange", "red","pink","darkblue","blue", "brown", "green",  "grey"))+ 
#  theme_set(theme_bw(base_size=12))

#ggplot (count_i, aes(x=PC5,y=PC6,col=stage,shape=sex)) + #this does not give variance percentage.
#  geom_point(size=3,alpha=0.8)+
#  scale_color_manual(values = c("orange", "red","pink","darkblue","blue", "brown", "green",  "grey"))+ 
#  theme_set(theme_bw(base_size=12)) +
#  theme(legend.justification=c(1,0), legend.position=c(0.95,0.05))

#pca1 <- PCA(pcaData)
#par(mar=c(5,5,4,3))
#plot.PCA(pca1, axes=c(1, 2), label="none", col.ind=c("orange", "orange", "orange", "orange","orange","orange","red","red","red","red","red","red","pink","pink", "pink", "pink", "pink","pink","pink","pink","pink","darkblue","darkblue", "darkblue", "darkblue","darkblue","darkblue","blue","blue","blue","blue","blue", "brown", "brown", "brown", "brown", "brown", "brown", "brown", "brown", "brown","brown", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"), title="PCA 1 vs PCA 2", cex=1.5)
#plotellipses(plot.PCA)
#plot.PCA(pca1, axes=c(3, 4), col.ind=c("red","red","red","blue","blue","blue","red","red","red","blue","blue","blue","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","blue","blue","blue","red","red","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue"), title="PCA 3 vs PCA 4", cex=0.7)
#plot.PCA(pca1, axes=c(4, 5), col.ind=c("red","red","red","blue","blue","blue","red","red","red","blue","blue","blue","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","blue","blue","blue","red","red","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue","blue","red","red","red","red","red","blue","blue","blue","blue"), title="PCA 1 vs PCA 2", cex=0.4)
#plotellipses(pca1)

#y <- dgl
#colnames(y) <- paste(colnames(y), design$group, sep="\n")
#cols = c(col.M23, col.F23,col.M27,col.F27, col.M31, col.F31,col.M43,col.F43,col.R43,col.M46,col.F46,col.R46)#for tvedora
#cols = c(col.M23, col.F23,col.U23,col.M27,col.F27,col.U27, col.M31, col.F31,col.U31,col.M43,col.F43,col.M46,col.F46,col.MF46) #for argovie
#cols = c(col.M23, col.F23,col.M27,col.F27, col.M31, col.F31,col.M43,col.F43,col.M46,col.F46,col.SR46,col.FB,col.BM,col.BF,col.FL,col.ML,col.FO,col.MT) #for ammarnas
#pchs = c(18,5,18,5,18,5,18,5,16,18,5,16) #for tvedora
#pchs = c(18,5,18,5,18,5,18,5,18,5,1,16,1,16,1,16,1)
#plotMDS(y, pch=pchs[design$group],col=cols[design$group], cex=0.6, main="Ammarnas MDS plot",cex.main=0.8, cex.lab=0.5,lty=1.5, lwd=2)
#legend('bottomright', inset=0.02, legend=levels(design$group), pch = pchs, col=cols )
#plotMDS(y, cex=1, col=cols[design$group], main=paste(sub_analyse,"MDS plot"))
#plotMDS(y, cex=1, pch=as.numeric(y$samples$group), col=as.numeric(y$samples$group), main=paste(sub_analyse,"MDS plot"))
#legend('bottomright', legend=levels(design$group), pch = as.numeric(y$samples$group), col=as.numeric(y$samples$group))
#dev.off()
