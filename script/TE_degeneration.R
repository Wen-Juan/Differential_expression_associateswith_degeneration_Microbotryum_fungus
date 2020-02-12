#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
install.packages("picante")
library(picante)

##analyzing data with randomizing non-DE genes between a1 and a2.
TE_homolog <- read.table('~/input/TE_degeneration/Mvsl_a1a2_te_nonDE_compart.txt', header = T)
str(TE_homolog)

geneok <- data.frame(TE_homolog$genea10k, TE_homolog$genea20k)
up5k <- data.frame(TE_homolog$up5k, TE_homolog$a2up5k)
up10k <- data.frame(TE_homolog$upk10, TE_homolog$a2upk10)
up15k <- data.frame(TE_homolog$upk15, TE_homolog$a2upk15)
up20k <- data.frame(TE_homolog$upk20, TE_homolog$a2upk20)
down5k <- data.frame(TE_homolog$down5k, TE_homolog$a2down5k)
down10k <- data.frame(TE_homolog$down10k, TE_homolog$a2down10k)
down15k <- data.frame(TE_homolog$down15k, TE_homolog$a2down15k)
down20k <- data.frame(TE_homolog$down20k, TE_homolog$a2down20k)

str(up5k)

geneok_rand <- randomizeMatrix(geneok,null.model = "richness",iterations = 1000)
up5k_rand <- randomizeMatrix(up5k,null.model = "richness",iterations = 1000)
head(up5k_rand)
up10k_rand <- randomizeMatrix(up10k,null.model = "richness",iterations = 1000)
up15k_rand <- randomizeMatrix(up15k,null.model = "richness",iterations = 1000)
up20k_rand <- randomizeMatrix(up20k,null.model = "richness",iterations = 1000)
down5k_rand <- randomizeMatrix(down5k,null.model = "richness",iterations = 1000)
down10k_rand <- randomizeMatrix(down10k,null.model = "richness",iterations = 1000)
down15k_rand <- randomizeMatrix(down15k,null.model = "richness",iterations = 1000)
down20k_rand <- randomizeMatrix(down20k,null.model = "richness",iterations = 1000)

total <- cbind(geneok_rand,up5k_rand, up10k_rand, up15k_rand, up20k_rand, down5k_rand, down10k_rand,down15k_rand,down20k_rand)
head(total)

total_data <- cbind(TE_homolog, total)
write.table(total_data,file = "~/input/TE_degeneration/Mvsl_a1a2_te_exp_compart_low-high.txt",quote=F, row.names=T, sep='\t')

TE_homolog_mod <- read.table('~/input/TE_degeneration/Mvsl_a1a2_te_exp_compart_low-high.txt', header = T)
str(TE_homolog_mod)

##combine internal gene and upstream up to 10kb and TE insertions
pdf("/output/figures/Mvsl_TE_exp_correlation_geneandupstream10k_pooled_mod.pdf", width=8, height=8)
p_a <- ggplot(TE_homolog_mod, aes(x=genetoup10kdiff, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(size =2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  labs(x='TE insertion difference in genes', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_correlation_geneanddownstream10k_pooled_mod.pdf", width=8, height=8)
p_b <- ggplot(TE_homolog_mod, aes(x=genetodown10kdiff, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(size =2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,10) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  labs(x='TE insertion difference in genes', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_upanddownstream_new_2figures_bw.pdf", width=12, height=8)
par(mar=c(6,6,6,6))
plot_grid(p_a, p_b,labels=c('A','B'))
dev.off()

####individual distance interval.
pdf("~/output/figures/Mvsl_TE_exp_correlation_gene0k_pooled_mod.pdf", width=8, height=8)
p_a <- ggplot(TE_homolog_mod, aes(x=genediff, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(size =2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  labs(x='TE insertion difference in genes', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_correlation_up0-5k_pooled_mod.pdf", width=8, height=8)
p_b <- ggplot(TE_homolog_mod, aes(x=upk5diff, y=abs, color=DE2,shape=DE2)) +
  scale_shape_manual(values=c(16,1),guide=FALSE) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(size =2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  labs(x='TE insertion difference at upstream 0-5kb', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


pdf("~/output/figures/Mvsl_TE_exp_correlation_up5-10k_pooled_mod.pdf", width=8, height=8)
p_c <- ggplot(TE_homolog_mod, aes(x=upk10diff, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1), labels=c("DE","Non-DE")) +
  scale_color_manual(values = c("black","dark grey"),guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  theme(legend.position = c(0.25, 0.75)) +
  labs(x='TE insertion difference at upstream 5-10kb', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_upstream_new_3figures_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(p_a, p_b, p_c,labels=c('A','B','C'))
dev.off()


pdf("~/output/figures/Mvsl_TE_exp_correlation_up10-15k_pooled_mod.pdf", width=8, height=8)
p1 <- ggplot(TE_homolog_mod, aes(x=upk15diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  theme_bw() + 
  theme(legend.position = c(0.25, 0.75)) +
  labs(x='TE insertion difference at upstream 10-15kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("~/output/figures/Mvsl_TE_exp_correlation_up15-20k_pooled_mod.pdf", width=8, height=8)
p2 <-ggplot(TE_homolog_mod, aes(x=upk20diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  theme_bw() + 
  labs(x='TE insertion difference at upstream 15-20kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_correlation_down0-5k_pooled_mod.pdf", width=8, height=8)
p3 <-ggplot(TE_homolog_mod, aes(x=down5kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  theme_bw() + 
  ylim(0,13) + xlim(-13,13) +
  labs(x='TE insertion difference at downstream 0-5kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_correlation_down5-10k_pooled_mod.pdf", width=8, height=8)
p4 <- ggplot(TE_homolog_mod, aes(x=down10kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  labs(x='TE insertion difference at downstream 5-10kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

TE_homolog_mod_DE <- subset(TE_homolog_mod, TE_homolog_mod$DE2 == "DE")
pdf("~/output/figures/Mvsl_TE_exp_correlation_down10-15k_pooled_mod.pdf", width=8, height=8)
p5 <-ggplot(TE_homolog_mod_DE, aes(x=upk5diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  theme_bw() + 
  ylim(0,13) + xlim(-13,13) +
  labs(x='TE insertion difference at downstream 10-15kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down10-20k_pooled_mod.pdf", width=8, height=8)
p6 <-ggplot(TE_homolog_mod, aes(x=down20kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  labs(x='TE insertion difference at downstream 15-20kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("~/output/figures/Mvsl_TE_restupdownstream_6figures_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(p1,p2,p3,p4,p5,p6,labels=c('A','B','C','D','E','F'))
dev.off()

### pool a1 and a2 genes, for difference of TE insertion sites.
TE_homolog_mod_rmcentro <- subset(TE_homolog_mod, TE_homolog_mod$youngold != "Centro")
TE_homolog_mod_rmcentro_oldstrata <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "OldStrata")
TE_homolog_mod_rmcentro_youngstrata <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "ColorStrata")

pdf("~/output/figures/Mvsl_TE_exp_correlation_innergenes_pooled_mod.pdf", width=8, height=8)
p_a <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=genediff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"), guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,2) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.75)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='TE insertions in genes') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

wilcox.test(TE_homolog_mod_rmcentro_oldstrata$genediff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$genediff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4420, p-value = 0.9037
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk5diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk5diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 3654.5, p-value = 0.04011
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$down5kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$down5kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4405.5, p-value = 0.9711
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk10diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk10diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 3991.5, p-value = 0.2665
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$down10kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$down10kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4138, p-value = 0.4869
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk15diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk15diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4733, p-value = 0.3395
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$down15kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$down15kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4405, p-value = 0.9729
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$down20kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$down20kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4335.5, p-value = 0.8776
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk20diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk20diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4252.5, p-value = 0.6969

wilcox.test(TE_homolog_mod_rmcentro_youngstrata$genediff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$genediff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 45, p-value = NA
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk5diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk5diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 49, p-value = 0.7588
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$down5kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$down5kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 18, p-value = 0.04105
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk10diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk10diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 46.5, p-value = 0.9001
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$down10kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$down10kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 41.5, p-value = 0.8368
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk15diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk15diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 46.5, p-value = 0.942
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$down15kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$down15kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 45, p-value = 1
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk20diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk20diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 32, p-value = 0.3909
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$down20kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$down20kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 29.5, p-value = 0.3085

TE_homolog_mod_rmcentro_auto <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "Auto")
TE_homolog_mod_rmcentro_PAR <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "bPAR")
wilcox.test(TE_homolog_mod_rmcentro_auto$genediff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$genediff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#WW = 1959100, p-value = 0.3705
wilcox.test(TE_homolog_mod_rmcentro_auto$upk5diff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$upk5diff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1949100, p-value = 0.7441
wilcox.test(TE_homolog_mod_rmcentro_auto$down5kdiff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$down5kdiff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1973900, p-value = 0.01463
wilcox.test(TE_homolog_mod_rmcentro_auto$upk10diff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$upk10diff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1966200, p-value = 0.1011
wilcox.test(TE_homolog_mod_rmcentro_auto$down10kdiff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$down10kdiff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1969400, p-value = 0.05582
wilcox.test(TE_homolog_mod_rmcentro_auto$upk15diff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$upk15diff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1947800, p-value = 0.6447
wilcox.test(TE_homolog_mod_rmcentro_auto$down15kdiff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$down15kdiff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#= 1946300, p-value = 0.5654
wilcox.test(TE_homolog_mod_rmcentro_auto$down20kdiff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$down20kdiff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1948200, p-value = 0.7238
wilcox.test(TE_homolog_mod_rmcentro_auto$upk20diff[TE_homolog_mod_rmcentro_auto$DE2=='DE'],TE_homolog_mod_rmcentro_auto$upk20diff[TE_homolog_mod_rmcentro_auto$DE2=='NON'],exact = FALSE) 
#W = 1957500, p-value = 0.5354

wilcox.test(TE_homolog_mod_rmcentro_PAR$genediff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$genediff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 612, p-value = NA
wilcox.test(TE_homolog_mod_rmcentro_PAR$upk5diff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$upk5diff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W =  600, p-value = 0.6406
wilcox.test(TE_homolog_mod_rmcentro_PAR$down5kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$down5kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 612, p-value = NA
wilcox.test(TE_homolog_mod_rmcentro_PAR$upk10diff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$upk10diff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 618, p-value = 0.7532
wilcox.test(TE_homolog_mod_rmcentro_PAR$down10kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$down10kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 663, p-value = 0.003891
wilcox.test(TE_homolog_mod_rmcentro_PAR$upk15diff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$upk15diff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 600, p-value = 0.6405
wilcox.test(TE_homolog_mod_rmcentro_PAR$down15kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$down15kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 612, p-value = NA
wilcox.test(TE_homolog_mod_rmcentro_PAR$down20kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$down20kdiff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 714, p-value = 3.772e-05
wilcox.test(TE_homolog_mod_rmcentro_PAR$upk20diff[TE_homolog_mod_rmcentro_PAR$DE2=='DE'],TE_homolog_mod_rmcentro_PAR$upk20diff[TE_homolog_mod_rmcentro_PAR$DE2=='NON'],exact = FALSE) 
#W = 600, p-value = 0.6405

pdf("~/output/figures/Mvsl_TE_exp_cor_upstream0-5k_pooled_mod.pdf", width=8, height=8)
p_b <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=upk5diff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"),guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 0-5kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_upstream5-10k_pooled_mod.pdf", width=8, height=8)
p_c <-ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=upk10diff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"),guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 5-10kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_upstream10-15k_pooled_mod.pdf", width=8, height=8)
p_d <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=upk15diff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"), guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 10-15kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_upstream15-20k_pooled_mod.pdf", width=8, height=8)
p_e <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=upk20diff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"),guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 15-20kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_downstream0-5k_pooled_mod.pdf", width=8, height=8)
p_f <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=down5kdiff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"), guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Downstream 0-5kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_downstream5-10k_pooled_mod.pdf", width=8, height=8)
p_g <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=down10kdiff, fill=DE2)) +
  scale_fill_manual(values = c("white","dark grey"),guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Downstream 5-10kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_downstream10-15k_pooled_mod.pdf", width=8, height=8)
p_h <-ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=down15kdiff, fill=DE2)) +
  scale_fill_manual(values =  c("white","dark grey"),guide=FALSE) +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  theme_bw() +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Downstream 10-15kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("~/output/figures/Mvsl_TE_exp_cor_downstream15-20k_pooled_mod.pdf", width=8, height=8)
p_i <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=down20kdiff, fill=DE2)) +
  scale_fill_manual(values = c("white","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.75) +
  ylim(-10,10) +  
  theme_bw() + 
  theme(legend.position = c(0.2, 0.78)) +
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Downstream 15-20kb') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/MvslMvld_dn_ds_combine_new_3figures_bw.pdf", width=12, height=8)
par(mar=c(10,10,10,10))
plot_grid(p_e,p_d, p_c, p_b, p_a, p_f, p_g, p_h,p_i,labels=c('A','B','C','D','E','F','G','H','I'))
dev.off()