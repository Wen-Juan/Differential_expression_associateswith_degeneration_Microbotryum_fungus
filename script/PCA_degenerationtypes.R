# plotting PCA 
#install.packages('ggfortify')
#library('ggfortify')

#y <- dgl
#count_norm <- data.frame(y$counts)
#count_norm1 <- data.frame(t(count_norm[,1:61]))
#str(count_norm1)

#pcaData <- count_norm1
#pca <- prcomp(pcaData, scale. = TRUE)
#summary(pca)

#rownames(count_norm1) <- c("23F1", "23F2", "23F3", "23M1","23M3","23M2","27F1","27F2","27F3", "27M1","27M2","27M3","31F1","31F2", "31F3", "31F4", "31M1","31M2","31M3","31M4","31M5","43F1","43F2", "43F3", "43M1","43M2","43M3","46F1","46F2","46M1","46M2","46M3", "FB1", "FB2", "FB3", "FB4","FB5", "MB1", "MB2", "MB3", "MB4", "MB5", "FO1", "FO2", "FO3", "FO4", "FO5", "MT1", "MT2", "MT3", "MT4", "MT5", "FL1", "FL2", "FL3", "FL4", "FL5", "ML1", "ML2", "ML3", "ML4")
#str(count_norm1)
#count_norm1$stage <- c("23", "23", "23", "23","23","23","27","27","27", "27","27","27","31","31", "31", "31", "31","31","31","31","31","43","43", "43", "43","43","43","46","46","46","46","46","Brain","Brain","Brain","Brain","Brain","Brain","Brain","Brain","Brain","Brain","Gonad","Gonad","Gonad","Gonad","Gonad","Gonad","Gonad","Gonad","Gonad","Gonad","Liver","Liver","Liver","Liver","Liver","Liver","Liver","Liver","Liver")
#count_norm1$sex <- c("female","female","female","male","male","male","female","female","female","male","male","male","female","female","female","female","male","male","male","male","male","female","female","female","male","male","male","female","female","male","male","male","female","female","female","female","female","male","male","male","male","male","female","female","female","female","female","male","male","male","male","male","female","female","female","female","female","male","male","male","male")

#count_i <- data.frame(pca$x, stage=count_norm1$stage, sex=count_norm1$sex)
#str(count_i)

#pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/alltissues_pc1pc2.pdf")
#autoplot(pca, data = count_i, x=1, y=2, colour = 'stage', shape='sex', size=4, ylim=c(-0.3,0.5), xlim=c(-0.3,0.3))
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
