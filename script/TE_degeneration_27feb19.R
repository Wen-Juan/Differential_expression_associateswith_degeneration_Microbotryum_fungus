#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
install.packages("picante")
library(picante)

##analyzing data with randomizing non-DE genes between a1 and a2, on 12feb.2019
TE_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/12feb2019/Mvsl_a1a2_te_nonDE_compart.txt', header = T)
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
write.table(total_data,file = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/12feb2019/Mvsl_a1a2_te_nonDE_compart_randomdized3.txt",quote=F, row.names=T, sep='\t')

###12feb2019
TE_homolog_mod <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/12feb2019/Mvsl_a1a2_te_all_compart_mod.txt', header = T)
str(TE_homolog_mod)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_gene0k_pooled_mod.pdf", width=8, height=8)
p_a <- ggplot(TE_homolog_mod, aes(x=gene0kdiff, y=abs, color=DE2, shape=DE2)) +
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

y_m1 <- lm(abs ~ DE2/gene0kdiff-1, data = TE_homolog_mod)
summary(y_m1)
#####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE             1.887440   0.018095 104.308  < 2e-16 ***
  DE2NON            0.214049   0.004949  43.251  < 2e-16 ***
  DE2DE:gene0kdiff  0.108022   0.033461   3.228  0.00125 ** 
  DE2NON:gene0kdiff 0.004697   0.020738   0.226  0.82084    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4414 on 8545 degrees of freedom
Multiple R-squared:  0.5989,	Adjusted R-squared:  0.5987 
#####

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up0-5k_pooled_mod.pdf", width=8, height=8)
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

y_m2 <- lm(abs ~ DE2/upk5diff-1, data = TE_homolog_mod)
summary(y_m2)
#########
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.887906   0.018085 104.390  < 2e-16 ***
  DE2NON           0.214075   0.004946  43.281  < 2e-16 ***
  DE2DE:upk5diff   0.123460   0.028126   4.389 1.15e-05 ***
  DE2NON:upk5diff -0.012632   0.012452  -1.014     0.31    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4411 on 8545 degrees of freedom
Multiple R-squared:  0.5993,	Adjusted R-squared:  0.5992  
##########

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up5-10k_pooled_mod.pdf", width=8, height=8)
p_c <- ggplot(TE_homolog_mod, aes(x=upk10diff, y=abs, color=DE2, shape=DE2)) +
  scale_shape_manual(values=c(16,1), labels=c("DE","Non-DE"), name = "Bias direction") +
  scale_color_manual(values = c("black","dark grey"),guide=FALSE) +
  geom_point(size = 2.5) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  theme(legend.position = c(0.25, 0.75)) +
  labs(x='TE insertion difference at upstream 5-10kb', y='Gene expression ratio (|Log2(A1/A2)|)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_upstream_new_3figures_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(p_a, p_b, p_c,labels=c('A','B','C'))
dev.off()

y_m3 <- lm(abs ~ DE2/upk10diff-1, data = TE_homolog_mod)
summary(y_m3)
###########
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.890732   0.018108 104.416  < 2e-16 ***
  DE2NON           0.214028   0.004947  43.264  < 2e-16 ***
  DE2DE:upk10diff  0.114492   0.027340   4.188 2.85e-05 ***
  DE2NON:upk10diff 0.009655   0.016995   0.568     0.57    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4412 on 8545 degrees of freedom
Multiple R-squared:  0.5992,	Adjusted R-squared:  0.599 
F-statistic:  3194 on 4 and 8545 DF,  p-value: < 2.2e-16
#############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up10-15k_pooled_mod.pdf", width=8, height=8)
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

y_m4 <- lm(abs ~ DE2/upk15diff-1, data = TE_homolog_mod)
summary(y_m4)
##############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE             1.886288   0.018110 104.157   <2e-16 ***
  DE2NON            0.214054   0.004951  43.233   <2e-16 ***
  DE2DE:upk15diff  -0.042650   0.028108  -1.517    0.129    
DE2NON:upk15diff  0.012156   0.017214   0.706    0.480    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8545 degrees of freedom
Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5983 
F-statistic:  3185 on 4 and 8545 DF,  p-value: < 2.2e-16
#############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up15-20k_pooled_mod.pdf", width=8, height=8)
p2 <-ggplot(TE_homolog_mod, aes(x=upk20diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  theme_bw() + 
  labs(x='TE insertion difference at upstream 15-20kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m5 <- lm(abs ~ DE2/upk20diff-1, data = TE_homolog_mod)
summary(y_m5)
########################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.8872491  0.0181044 104.243   <2e-16 ***
  DE2NON           0.2140560  0.0049517  43.229   <2e-16 ***
  DE2DE:upk20diff  0.0205518  0.0164466   1.250    0.211    
DE2NON:upk20diff 0.0002931  0.0124611   0.024    0.981    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8545 degrees of freedom
Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5983 
F-statistic:  3184 on 4 and 8545 DF,  p-value: < 2.2e-16
#######################

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down0-5k_pooled_mod.pdf", width=8, height=8)
p3 <-ggplot(TE_homolog_mod, aes(x=down5kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  theme_bw() + 
  ylim(0,13) + xlim(-13,13) +
  labs(x='TE insertion difference at downstream 0-5kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m6 <- lm(abs ~ DE2/down5kdiff-1, data = TE_homolog_mod)
summary(y_m6)
#############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE              1.886071   0.018088 104.274  < 2e-16 ***
  DE2NON             0.214057   0.004947  43.273  < 2e-16 ***
  DE2DE:down5kdiff  -0.099744   0.023384  -4.266 2.02e-05 ***
  DE2NON:down5kdiff  0.006780   0.016993   0.399     0.69    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4412 on 8545 degrees of freedom
Multiple R-squared:  0.5992,	Adjusted R-squared:  0.5991 
F-statistic:  3194 on 4 and 8545 DF,  p-value: < 2.2e-16
#############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down5-10k_pooled_mod.pdf", width=8, height=8)
p4 <- ggplot(TE_homolog_mod, aes(x=down10kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-13,13) +
  theme_bw() + 
  labs(x='TE insertion difference at downstream 5-10kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m7 <- lm(abs ~ DE2/down10kdiff-1, data = TE_homolog_mod)
summary(y_m7)
##########
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE               1.88237    0.01816 103.656  < 2e-16 ***
  DE2NON              0.21403    0.00495  43.239  < 2e-16 ***
  DE2DE:down10kdiff  -0.06216    0.02020  -3.078  0.00209 ** 
  DE2NON:down10kdiff -0.00439    0.01551  -0.283  0.77717    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4414 on 8545 degrees of freedom
Multiple R-squared:  0.5988,	Adjusted R-squared:  0.5986 
F-statistic:  3189 on 4 and 8545 DF,  p-value: < 2.2e-16
#############

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down10-15k_pooled_mod.pdf", width=8, height=8)
p5 <-ggplot(TE_homolog_mod, aes(x=down15kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("black","dark grey"), guide = FALSE) +
  geom_point(alpha=0.8) + geom_smooth(method = lm) +
  theme_bw() + 
  ylim(0,13) + xlim(-13,13) +
  labs(x='TE insertion difference at downstream 10-15kb', y='Gene expression ratio |Log2(A1/A2)|') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y_m8 <- lm(abs ~ DE2/down15kdiff-1, data = TE_homolog_mod)
summary(y_m8)
###################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE               1.887235   0.018122 104.141   <2e-16 ***
  DE2NON              0.214058   0.004952  43.227   <2e-16 ***
  DE2DE:down15kdiff   0.005254   0.025521   0.206    0.837    
DE2NON:down15kdiff -0.002450   0.015357  -0.160    0.873    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8545 degrees of freedom
Multiple R-squared:  0.5984,	Adjusted R-squared:  0.5982 
F-statistic:  3183 on 4 and 8545 DF,  p-value: < 2.2e-16 
###################

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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_restupdownstream_6figures_bw.pdf", width=12, height=8)
par(mar=c(10,10,6,6))
plot_grid(p1,p2,p3,p4,p5,p6,labels=c('A','B','C','D','E','F'))
dev.off()


y_m9 <- lm(abs ~ DE2/down20kdiff-1, data = TE_homolog_mod)
summary(y_m9)

##############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE               1.887160   0.018106 104.228   <2e-16 ***
  DE2NON              0.214053   0.004952  43.226   <2e-16 ***
  DE2DE:down20kdiff   0.012485   0.025584   0.488    0.626    
DE2NON:down20kdiff -0.003213   0.015404  -0.209    0.835   
##############

### pool a1 and a2 genes, for difference of TE insertion sites, analysis on 26 Feb.
TE_homolog_mod_rmcentro <- subset(TE_homolog_mod, TE_homolog_mod$youngold != "Centro")
TE_homolog_mod_rmcentro_oldstrata <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "OldStrata")
TE_homolog_mod_rmcentro_youngstrata <- subset(TE_homolog_mod, TE_homolog_mod$youngold == "ColorStrata")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_innergenes_pooled_mod.pdf", width=8, height=8)
p_a <- ggplot(TE_homolog_mod_rmcentro, aes(x=youngold, y=gene0kdiff, fill=DE2)) +
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

wilcox.test(TE_homolog_mod_rmcentro_oldstrata$gene0kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$gene0kdiff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4538, p-value = 0.5221
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk5diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk5diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4436, p-value = 0.9036
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk10diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk10diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4196, p-value = 0.5576
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk15diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk15diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4479, p-value = 0.8083
wilcox.test(TE_homolog_mod_rmcentro_oldstrata$upk20diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='DE'],TE_homolog_mod_rmcentro_oldstrata$upk20diff[TE_homolog_mod_rmcentro_oldstrata$DE2=='NON'],exact = FALSE) 
#W = 4398.5, p-value = 0.9866

wilcox.test(TE_homolog_mod_rmcentro_youngstrata$gene0kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$gene0kdiff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 45, p-value = NA
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk5diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk5diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 45.5, p-value = 1
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk10diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk10diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 46.5, p-value = 0.9001
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk15diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk15diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 49.5, p-value = 0.7707
wilcox.test(TE_homolog_mod_rmcentro_youngstrata$upk20diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='DE'],TE_homolog_mod_rmcentro_youngstrata$upk20diff[TE_homolog_mod_rmcentro_youngstrata$DE2=='NON'],exact = FALSE) 
#W = 38, p-value = 0.6555

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_upstream0-5k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_upstream5-10k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_upstream10-15k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_upstream15-20k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_downstream0-5k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_downstream5-10k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_downstream10-15k_pooled_mod.pdf", width=8, height=8)
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_cor_downstream15-20k_pooled_mod.pdf", width=8, height=8)
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

#loading data on 01Feb.2019.
##below are genes with single copy homologs.
TE_number_touse_singlecopy_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/TE_number_touse_01feb19.txt', header = T)
str(TE_number_touse_singlecopy_homolog)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_prop_all_intervals_01feb19.pdf", width=8, height=8)
ggplot(TE_number_touse_singlecopy_homolog, aes(x=location,y=Te_prop, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                
  theme_bw() + 
  theme(legend.position = c(0.2, 0.85)) +
  scale_x_discrete(labels=c("up:15-20k", "10-15k","5-10k","0-5k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval windows") +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/geneswithte_prop_all_intervals_01feb19.pdf", width=8, height=8)
ggplot(TE_number_touse_singlecopy_homolog, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of genes with TE insertions', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

### loading dataset of individual gene with TE insertions, 03Feb.2019
TE_singlegene <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/Mvsl_a1a2_te_exp_compart.txt', header = T)
str(TE_singlegene)

TE_singlegene_rmcentro <- subset(TE_singlegene, TE_singlegene$youngold!="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genediff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=genediff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-2,2) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Gene difference in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_geneup5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down5kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown10kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down10kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown15kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down15kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_genedown20kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=down20kdiff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-10,10) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up10kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 5-10k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up15kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 10-15k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_up20kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=youngold, y=upk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-5,5) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='Upstream 15-20k in TE insertions (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_gene0k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=genediff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for homologs introns (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_gene0k_pooled.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=genediff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for homologs introns (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_gene0k_pooled.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=genediff, y=abs)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for homologs introns (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up5k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk5diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-15,15) +
  labs(x='TE difference for upstream 0-5k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down5k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down5kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for downstream 0-5k (A1-A2)', y='Absolute value of Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down5k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down5kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 0-5k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up10k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk10diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 5-10k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up10k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk10diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for upstream 5-10k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down10k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down10kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 5-10k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down10k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down10kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for downstream 5-10k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down15k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down15kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 10-105k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down15k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down15kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for downstream 10-15k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down20k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down20kdiff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for downstream 15-20k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_down20k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=down20kdiff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for downstream 15-20k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up15k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk15diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 10-15k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up15k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk15diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for upstream 10-15k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up20k.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk20diff, y=logFC.A1.A2, color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(-10,10) + xlim(-12,12) +
  labs(x='TE difference for upstream 15-20k (A1-A2)', y='Gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_TE_exp_correlation_up20k_depool.pdf", width=8, height=8)
ggplot(TE_singlegene_rmcentro, aes(x=upk20diff, y=abs, color=DE2)) +
  scale_color_manual(values = c("firebrick3","dark grey"),labels=c("DE","Non-DE"), name = "Bias direction") +
  geom_point() + geom_smooth(method = lm) +
  ylim(0,13) + xlim(-10,10) +
  labs(x='TE difference for upstream 15-20k (A1-A2)', y='Absolute value of gene expression ratio Log2(A1/A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###stats.12feb2019
y1 <- lm(abs ~ DE2/genediff-1, data = TE_singlegene_rmcentro)
summary(y1)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.889250   0.018118 104.275   <2e-16 ***
  DE2NON           0.213934   0.004952  43.202   <2e-16 ***
  DE2DE:genediff  -0.057360   0.033572  -1.709   0.0876 .  
DE2NON:genediff  0.018541   0.020750   0.894   0.3716    
---
Residual standard error: 0.4415 on 8544 degrees of freedom
Multiple R-squared:  0.5987,	Adjusted R-squared:  0.5985 
F-statistic:  3186 on 4 and 8544 DF,  p-value: < 2.2e-16
##

y1a <- lm(abs ~ DE/genediff-1, data = TE_singlegene_rmcentro)
summary(y1a)

###################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown           1.768205   0.030942  57.146  < 2e-16 ***
  DENON            0.213934   0.004943  43.280  < 2e-16 ***
  DEUp             1.951272   0.022291  87.538  < 2e-16 ***
  DEDown:genediff  0.213587   0.089989   2.373  0.01764 *  
  DENON:genediff   0.018541   0.020713   0.895  0.37072    
DEUp:genediff   -0.101458   0.036109  -2.810  0.00497 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4407 on 8542 degrees of freedom
Multiple R-squared:  0.6002,	Adjusted R-squared:  0.5999 
F-statistic:  2137 on 6 and 8542 DF,  p-value: < 2.2e-16
###################################

y2 <- lm(abs ~ DE2/upk5diff-1, data = TE_singlegene_rmcentro)
summary(y2)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.888035   0.018145 104.054   <2e-16 ***
  DE2NON           0.214060   0.004951  43.236   <2e-16 ***
  DE2DE:upk5diff  -0.015253   0.028253  -0.540    0.589    
DE2NON:upk5diff  0.002611   0.012464   0.209    0.834    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8544 degrees of freedom
Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5983 
F-statistic:  3184 on 4 and 8544 DF,  p-value: < 2.2e-16
##
y2a <- lm(abs ~ DE/upk5diff-1, data = TE_singlegene_rmcentro)
summary(y2a)

##############################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown           1.775244   0.031003  57.260  < 2e-16 ***
  DENON            0.214060   0.004940  43.330  < 2e-16 ***
  DEUp             1.946515   0.022303  87.275  < 2e-16 ***
  DEDown:upk5diff  0.083467   0.037603   2.220 0.026463 *  
  DENON:upk5diff   0.002611   0.012437   0.210 0.833752    
DEUp:upk5diff   -0.150306   0.042635  -3.525 0.000425 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4406 on 8542 degrees of freedom
Multiple R-squared:  0.6003,	Adjusted R-squared:  0.6001 
F-statistic:  2138 on 6 and 8542 DF,  p-value: < 2.2e-16
######################

y2 <- lm(abs ~ DE2/upk10diff-1, data = TE_singlegene_rmcentro)
summary(y2)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE             1.888352   0.018116 104.235   <2e-16 ***
  DE2NON            0.214010   0.004951  43.227   <2e-16 ***
  DE2DE:upk10diff  -0.032902   0.027814  -1.183    0.237    
DE2NON:upk10diff  0.011818   0.014221   0.831    0.406    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4415 on 8544 degrees of freedom
Multiple R-squared:  0.5986,	Adjusted R-squared:  0.5984 
F-statistic:  3185 on 4 and 8544 DF,  p-value: < 2.2e-16
##

y2b <- lm(abs ~ DE/upk10diff-1, data = TE_singlegene_rmcentro)
summary(y2b)
##################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown            1.775283   0.031005  57.258  < 2e-16 ***
  DENON             0.214010   0.004938  43.335  < 2e-16 ***
  DEUp              1.954177   0.022292  87.661  < 2e-16 ***
  DEDown:upk10diff  0.077660   0.037741   2.058   0.0397 *  
  DENON:upk10diff   0.011818   0.014185   0.833   0.4048    
DEUp:upk10diff   -0.182547   0.041105  -4.441 9.07e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4404 on 8542 degrees of freedom
Multiple R-squared:  0.6007,	Adjusted R-squared:  0.6004 
F-statistic:  2142 on 6 and 8542 DF,  p-value: < 2.2e-16
####################################

y3 <- lm(abs ~ DE2/upk15diff-1, data = TE_singlegene_rmcentro)
summary(y3)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE             1.8893486  0.0181079 104.338  < 2e-16 ***
  DE2NON            0.2140582  0.0049482  43.259  < 2e-16 ***
  DE2DE:upk15diff  -0.0920261  0.0283112  -3.251  0.00116 ** 
  DE2NON:upk15diff -0.0004103  0.0172040  -0.024  0.98097    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4413 on 8544 degrees of freedom
Multiple R-squared:  0.599,	Adjusted R-squared:  0.5988 
F-statistic:  3190 on 4 and 8544 DF,  p-value: < 2.2e-16
##

y3a <- lm(abs ~ DE/upk15diff-1, data = TE_singlegene_rmcentro)
summary(y3a)
################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown            1.7688178  0.0309307  57.186  < 2e-16 ***
  DENON             0.2140582  0.0049412  43.321  < 2e-16 ***
  DEUp              1.9507456  0.0222960  87.493  < 2e-16 ***
  DEDown:upk15diff -0.1514183  0.0409175  -3.701 0.000216 ***
  DENON:upk15diff  -0.0004103  0.0171794  -0.024 0.980945    
DEUp:upk15diff   -0.0430603  0.0391213  -1.101 0.271064    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4407 on 8542 degrees of freedom
Multiple R-squared:  0.6002,	Adjusted R-squared:  0.5999 
F-statistic:  2137 on 6 and 8542 DF,  p-value: < 2.2e-16
################################

y4 <- lm(abs ~ DE2/upk20diff-1, data = TE_singlegene_rmcentro)
summary(y4)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE            1.888465   0.018142 104.093   <2e-16 ***
  DE2NON           0.214060   0.004951  43.235   <2e-16 ***
  DE2DE:upk20diff  0.001909   0.016478   0.116    0.908    
DE2NON:upk20diff 0.002060   0.012459   0.165    0.869    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8544 degrees of freedom
Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5983 
F-statistic:  3184 on 4 and 8544 DF,  p-value: < 2.2e-16
##
y4a <- lm(abs ~ DE/upk20diff-1, data = TE_singlegene_rmcentro)
summary(y4a)
#################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown            1.766350   0.031047  56.893   <2e-16 ***
  DENON             0.214060   0.004944  43.294   <2e-16 ***
  DEUp              1.950642   0.022319  87.397   <2e-16 ***
  DEDown:upk20diff  0.053589   0.033249   1.612    0.107    
DENON:upk20diff   0.002060   0.012442   0.166    0.869    
DEUp:upk20diff   -0.013740   0.018940  -0.725    0.468    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.441 on 8542 degrees of freedom
Multiple R-squared:  0.5997,	Adjusted R-squared:  0.5994 
F-statistic:  2133 on 6 and 8542 DF,  p-value: < 2.2e-16
###############################

y5 <- lm(abs ~ DE2/down5kdiff-1, data = TE_singlegene_rmcentro)
summary(y5)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE              1.892999   0.018058 104.827  < 2e-16 ***
  DE2NON             0.214063   0.004933  43.397  < 2e-16 ***
  DE2DE:down5kdiff  -0.187756   0.023458  -8.004 1.37e-15 ***
  DE2NON:down5kdiff  0.005544   0.016945   0.327    0.744    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4399 on 8544 degrees of freedom
Multiple R-squared:  0.6015,	Adjusted R-squared:  0.6013 
F-statistic:  3224 on 4 and 8544 DF,  p-value: < 2.2e-16
##
y5a <- lm(abs ~ DE/down5kdiff-1, data = TE_singlegene_rmcentro)
summary(y5a)
####################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown             1.777670   0.030812  57.694   <2e-16 ***
  DENON              0.214063   0.004921  43.503   <2e-16 ***
  DEUp               1.951875   0.022206  87.898   <2e-16 ***
  DEDown:down5kdiff -0.298816   0.033571  -8.901   <2e-16 ***
  DENON:down5kdiff   0.005544   0.016904   0.328   0.7429    
DEUp:down5kdiff   -0.082573   0.032638  -2.530   0.0114 *  
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4388 on 8542 degrees of freedom
Multiple R-squared:  0.6035,	Adjusted R-squared:  0.6032 
F-statistic:  2167 on 6 and 8542 DF,  p-value: < 2.2e-16
###################################

y6 <- lm(abs ~ DE2/down10kdiff-1, data = TE_singlegene_rmcentro)
summary(y6)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE              1.8881392  0.0180741 104.467  < 2e-16 ***
  DE2NON             0.2140572  0.0049392  43.339  < 2e-16 ***
  DE2DE:down10kdiff  0.1291383  0.0201061   6.423 1.41e-10 ***
  DE2NON:down10kdiff 0.0001266  0.0154777   0.008    0.993    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4405 on 8544 degrees of freedom
Multiple R-squared:  0.6004,	Adjusted R-squared:  0.6002 
F-statistic:  3209 on 4 and 8544 DF,  p-value: < 2.2e-16
##

y6a <- lm(abs ~ DE/down10kdiff-1, data = TE_singlegene_rmcentro)
summary(y6a)
################################
Residuals:
  Min      1Q  Median      3Q     Max 
-2.1096 -0.1538 -0.0674  0.0821  9.7737 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown             1.7756899  0.0309953  57.289  < 2e-16 ***
  DENON              0.2140572  0.0049307  43.413  < 2e-16 ***
  DEUp               1.9380755  0.0223066  86.883  < 2e-16 ***
  DEDown:down10kdiff 0.0520091  0.0275472   1.888   0.0591 .  
DENON:down10kdiff  0.0001266  0.0154511   0.008   0.9935    
DEUp:down10kdiff   0.2022874  0.0295372   6.849 7.98e-12 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4397 on 8542 degrees of freedom
Multiple R-squared:  0.6019,	Adjusted R-squared:  0.6016 
F-statistic:  2152 on 6 and 8542 DF,  p-value: < 2.2e-16
#################################

y7 <- lm(abs ~ DE2/down15kdiff-1, data = TE_singlegene_rmcentro)
summary(y7)
##
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE              1.888557   0.018117 104.240   <2e-16 ***
  DE2NON             0.214058   0.004951  43.235   <2e-16 ***
  DE2DE:down15kdiff  0.010144   0.025536   0.397    0.691    
DE2NON:down15kdiff 0.003114   0.015354   0.203    0.839    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4416 on 8544 degrees of freedom
Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5983 
F-statistic:  3184 on 4 and 8544 DF,  p-value: < 2.2e-16
##

y7a <- lm(abs ~ DE/down15kdiff-1, data = TE_singlegene_rmcentro)
summary(y7a)
################################
Residuals:
  Min      1Q  Median      3Q     Max 
-1.3889 -0.1538 -0.0676  0.0821  9.7619 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown             1.770661   0.031002  57.115   <2e-16 ***
  DENON              0.214058   0.004945  43.286   <2e-16 ***
  DEUp               1.949903   0.022316  87.376   <2e-16 ***
  DEDown:down15kdiff 0.008920   0.043313   0.206    0.837    
DENON:down15kdiff  0.003114   0.015336   0.203    0.839    
DEUp:down15kdiff   0.003128   0.031600   0.099    0.921    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.441 on 8542 degrees of freedom
Multiple R-squared:  0.5995,	Adjusted R-squared:  0.5992 
F-statistic:  2131 on 6 and 8542 DF,  p-value: < 2.2e-16
################################

y8 <- lm(abs ~ DE2/down20kdiff-1, data = TE_singlegene_rmcentro)
summary(y8)

##
Residuals:
  Min      1Q  Median      3Q     Max 
-1.4505 -0.1541 -0.0674  0.0819  9.8235 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DE2DE               1.888263   0.018112 104.256   <2e-16 ***
  DE2NON              0.214069   0.004950  43.250   <2e-16 ***
  DE2DE:down20kdiff  -0.061644   0.025614  -2.407   0.0161 *  
  DE2NON:down20kdiff -0.003821   0.015397  -0.248   0.8040    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4414 on 8544 degrees of freedom
Multiple R-squared:  0.5988,	Adjusted R-squared:  0.5986 
F-statistic:  3187 on 4 and 8544 DF,  p-value: < 2.2e-16
##
y8a <- lm(abs ~ DE/down20kdiff-1, data = TE_singlegene_rmcentro)
summary(y8a)

#########################################
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
DEDown              1.769518   0.030948  57.177   <2e-16 ***
  DENON               0.214069   0.004944  43.302   <2e-16 ***
  DEUp                1.949975   0.022296  87.457   <2e-16 ***
  DEDown:down20kdiff -0.053598   0.035884  -1.494   0.1353    
DENON:down20kdiff  -0.003821   0.015378  -0.248   0.8038    
DEUp:down20kdiff   -0.072406   0.036488  -1.984   0.0472 *  
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4409 on 8542 degrees of freedom
Multiple R-squared:  0.5998,	Adjusted R-squared:  0.5995 
F-statistic:  2134 on 6 and 8542 DF,  p-value: < 2.2e-16
#########################################

###stats.01feb2019
y <- lm(logFC.A1.A2 ~ upk5diff*DE, data = TE_singlegene_rmcentro)
summary(y)
y1 <- lm(logFC.A1.A2 ~ upk5diff+DE, data = TE_singlegene_rmcentro)
anova(y,y1)
###
Model 1: logFC.A1.A2 ~ upk5diff * DE
Model 2: logFC.A1.A2 ~ upk5diff + DE
Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
1   8542 2022.3                                
2   8544 2025.0 -2   -2.7503 5.8085 0.003014 **
###
 
y2 <- lm(logFC.A1.A2 ~ upk10diff*DE, data = TE_singlegene_rmcentro)
summary(y2)
y3 <- lm(logFC.A1.A2 ~ upk10diff+DE, data = TE_singlegene_rmcentro)
anova(y2,y3) 
###

y6 <- lm(logFC.A1.A2 ~ down5kdiff*DE, data = TE_singlegene_rmcentro)
summary(y6)
y7 <- lm(logFC.A1.A2 ~ down5kdiff+DE, data = TE_singlegene_rmcentro)
anova(y6,y7) 

TE_singlegene_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/01feb19/Mvsl_a1a2_te_exp_compart_sep.txt', header = T)
str(TE_singlegene_sep)

TE_singlegene_sep_rmcentro <- subset(TE_singlegene_sep, TE_singlegene_sep$youngold != "Centro")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_geneup5kdiff_exp.pdf", width=8, height=8)
ggplot(TE_singlegene_sep_rmcentro, aes(x=youngold, y=genea10k, fill=DE2)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm(logFC.A1.A2 ~ genediff*DE, data = TE_singlegene_rmcentro)
summary(y)
y1 <- lm(logFC.A1.A2 ~ genediff+DE, data = TE_singlegene_rmcentro)
anova(y,y1)
###

##below are all genes
TE_number_touse <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/TE_number_touse.txt', header = T)
str(TE_number_touse)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_prop_all_intervals.pdf", width=8, height=8)
ggplot(TE_number_touse, aes(x=location,y=Te_prop, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.2) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/genes_withTEinsertions_intervals.pdf", width=8, height=8)
ggplot(TE_number_touse, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.2) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","gene","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#loading data on 24Jan.2019, check distribution of coding sequence length
a1cds_length <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/A1_cds_length.txt', header = F)
str(a1cds_length)
a2cds_length <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/A2_cds_length.txt', header = F)
str(a2cds_length)
length_diff <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/A1A2_homolog/homolog_length_variation.txt', header = F)
str(length_diff)

ggplot(data=a1cds_length, aes(a1cds_length$V2)) + 
  xlim(0,7500) +
  ylim(0,2000) +
  geom_histogram(binwidth=200)

ggplot(data=a2cds_length, aes(a2cds_length$V2)) + 
  xlim(0,7500) +
  geom_histogram(binwidth=200)

ggplot(length_diff, aes(length_diff$V8)) + 
  xlim(0,1000) +
  ylim(0,40) +
  geom_histogram(binwidth=25)

#load the corresponding data files. modify this on Jan.17.2019.
#summary of TE numbers
TE_number <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/TE_nr_17012019.txt', header = T)
str(TE_number)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_number_intervals_5k.pdf", width=8, height=8)
ggplot(TE_number, aes(x=location,y=TE_interval, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,11000) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Number of detected TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Prop_geneswith_TEinsertion_intervals_5k.pdf", width=8, height=8)
ggplot(TE_number, aes(x=location,y=prop_withte, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.32) +                    
  scale_x_discrete(labels=c("up:20-15k", "15-10k","10-5k","5-0k","down:0-5k", "5-10k","10-15k","15-20k")) + 
  labs(y='Proportion of total genes with TE insertions', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#load the corresponding data files, DE genes and TE intervals. modify this on Jan.17.2019.
TE_DEexp <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/Mvsl_TE_exp_allintervals.txt', header = T)
str(TE_DEexp)

TE_DEexp1 <- subset(TE_DEexp, TE_DEexp$youngold != "Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y <- lm (logFC.A1.A2~upk20diff*DE-1, data=TE_DEexp1)
summary(y)

###
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
upk5diff       -0.159073   0.041431  -3.840 0.000124 ***
  DEDown         -1.851788   0.031855 -58.132  < 2e-16 ***
  DENON          -0.004183   0.004791  -0.873 0.382622    
DEUp            1.819751   0.036815  49.430  < 2e-16 ***
  upk5diff:DENON  0.192189   0.042263   4.547 5.51e-06 ***
  upk5diff:DEUp   0.237904   0.051473   4.622 3.86e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4227 on 8348 degrees of freedom
Multiple R-squared:  0.4217,	Adjusted R-squared:  0.4213 
F-statistic:  1015 on 6 and 8348 DF,  p-value: < 2.2e-16
###

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk5diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='0-5k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5-10kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_5-10kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk10diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='5-10k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_10-15kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_10-15kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk15diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='10-15k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_15-20kup_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=upk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k upstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_15-20kdown_exp.pdf", width=8, height=8)
ggplot(TE_DEexp1, aes(x=youngold, y=downk20diff, fill=DE)) + 
  scale_fill_manual(values = c("firebrick3","light grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(-4,3) +  
  scale_x_discrete(labels=c("Autosome", "PAR", "Young strata","Old strata")) + 
  labs(x='Genomic compartment', y='15-20k downstream difference (A1-A2)') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###only working on data on MAT chromosome.
TE_DEexp_mat <- subset(TE_DEexp, TE_DEexp$youngold != "Auto")
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_0-5kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk5diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() +geom_smooth(method=lm) +
  ylim(-12,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 5k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_0-5kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk5diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-6,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 5k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_5-10kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk10diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 5-10k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_5-10kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk10diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 5-10k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_10-15kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk15diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 10-15k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_10-15kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk15diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 10-15k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_15-20kup_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=upk20diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Upstream 15-20k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_MATonly_15-20kdown_exp_scatter.pdf", width=8, height=8)
ggplot(TE_DEexp_mat, aes(x=logFC.A1.A2, y=downk20diff,color=DE)) +
  scale_color_manual(values = c("firebrick3","dark grey","dodgerblue3"),labels=c("A2-biased","Not-biased","A1-biased"), name="Bias direction") +
  geom_point() + geom_smooth(method=lm) +
  ylim(-13,6) +
  xlim(-13,13) +
  labs(x='Expression Log2(A1/A2)', y='Downstream 15-20k (A1-A2))') +
  theme(axis.title.x = element_text(size=12,colour = "black"),axis.title.y = element_text(size=12,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#load TE and exp. data separately, Jan.18.2018
UpTE_DEexp_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/UpstreamTE_Exp_intervals_sep.txt', header = T)
str(UpTE_DEexp_sep)

UpTE_DEexp_sep1 <- subset(UpTE_DEexp_sep, UpTE_DEexp_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_upstream_exp.pdf", width=8, height=8)
ggplot(UpTE_DEexp_sep1, aes(x=youngold, y=up5k, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +
  facet_grid(~Up) +
  scale_x_discrete(labels=c("Auto", "PAR", "Young","Old")) + 
  labs(x='Genomic compartment', y='Average TE insertion sites') +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

DownTE_DEexp_sep <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/16jan2019/to_makegrah/DownstreamTE_Exp_intervals_sep.txt', header = T)
str(DownTE_DEexp_sep)
DownTE_DEexp_sep1 <- subset(DownTE_DEexp_sep, DownTE_DEexp_sep$youngold !="Centro")

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TEinterval_downstream_exp.pdf", width=8, height=8)
ggplot(DownTE_DEexp_sep1, aes(x=youngold, y=downak5, fill=interaction(haploid,DE))) + 
  scale_fill_manual(values = c("firebrick2","firebrick4","light grey","dark grey","dodgerblue2","dodgerblue4"), labels=c("A2-biased at A1","A2-biased at A2","Not-biased at A1","Not-biased at A2","A1-biased at A1","A1-biased at A2"), name="Bias direction") +
  geom_boxplot(notch=FALSE,outlier.shape=NA,alpha=0.85) +
  ylim(0,6) +
  facet_grid(~Down) +
  scale_x_discrete(labels=c("Auto", "PAR", "Young","Old")) + 
  labs(x='Genomic compartment', y='TE insertion at 5k-upstream (A1-A2)') +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#load the corresponding data files. modify this on Jan.09.2019.
#summary of TE numbers
TE_number <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/TE_number.txt', header = T)
str(TE_number)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/TE_number_intervals.pdf", width=8, height=8)
ggplot(TE_number, aes(x=range, y=TE, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,4000) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Number of detected TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#prop. gene with TE insertions
Prop_gene <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/gene_TEinsertion_prop.txt', header = T)
str(Prop_gene)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Gene_proportion_intervals.pdf", width=8, height=8)
ggplot(Prop_gene, aes(x=range, y=proportion, fill=haploid)) + 
  scale_fill_manual(values = c("dodgerblue3","firebrick3"), labels=c("A1","A2"), name="Haploid") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(0,0.15) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of TE insertion sites', x="Interval window") +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#de genes on MAT 
exp_DE_MAT_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_mat_cor.txt', header = T)
str(exp_DE_MAT_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_diffMAT_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_MAT_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.25,0.4) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))

dev.off()


#DE autosome
exp_DE_auto_cor <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_auto_cor.txt', header = T)
str(exp_DE_auto_cor)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_autosomes_DEnonDEinterval_overlapwithTE_prop_genes_cor.pdf", width=8, height=8)
ggplot(exp_DE_auto_cor, aes(x=interval, y=diffprop, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.2,0.2) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

#################################not consider this at the moment################################################
### richmedium
DE_richmedium <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/DEexp_TE_richmedium.txt', header = T)
str(DE_richmedium)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDEinterval_overlapwithTE_propgenes_richmeidum_cor.pdf", width=8, height=8)
ggplot(DE_richmedium, aes(x=interval, y=diffpropcor, fill=bias)) + 
  scale_fill_manual(values = c("firebrick3","grey","dodgerblue3"), labels=c("Down","NON", "Up"), name="DE expression") +  
  geom_bar(stat="identity",position=position_dodge(),alpha=0.85,lwd=0.5) +
  ylim(-0.1,0.1) +                    
  scale_x_discrete(labels=c("up:20-10k", "2-10k","0-2k","down:0-2k", "2-10k","10-20k")) + 
  labs(y='Proportion of genes with TE insertion site') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
#############################################################################################

### 1,20k up and down stream TE insertion sites among genomic compartments.
DE_homolog <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_2k10kupdownstream.txt', header = T)
str(DE_homolog)

################### I think this figure is not correct, both catogory makes the boxplot look difference among DE and Non-DE genes.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
y<-ggplot(DE_homolog, aes(y=logFC.A1.A2, x=k2updiff, color=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), labels=c("DE","Non DE"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  facet_grid(~chr) +
  ylim(-4,4) +  
  xlim(-2,2) +
  labs(y='LogFC(A1/A2)') + 
  labs(x='Difference of TE number between homologs (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
y + coord_flip()
ddev.off()
############################################################################################################
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff2kdownstream_logfc_3genomicparts.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-0.5,0.5) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff2kupstream_logfc_3genomicparts,.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-0.5,0.5) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-2k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff8kupstream_logfc_3genomicparts", width=8, height=8)
  ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
    labs(x='Genomic compartment') + 
    labs(y='Difference of TE number [2-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()
  
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervalupdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10kdwonstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k10downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_8kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k8kintervaldowndiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2-10k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","grey", "dodgerblue2"), labels=c("Down", "NON", "Up"), name="DE expression") + 
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_2kupstream.pdf", width=8, height=8)
ggplot(DE_homolog, aes(y=k2updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-4,4) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [2k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

###
DE_homolog2 <- read.table('/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/input/TE_degeneration/Mvsl_DEnonDE_TEinsert_mod8.txt', header = T)
str(DE_homolog2)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20updiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=k20downdiff,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [0-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kupstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20updiffa1,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k upstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Haploidselection_and_dosagecompensation_in_Microbotryum/output/figures/Mvsl_DEnonDE_TEdiff_logfc_3genomicparts_10to20kdownstream.pdf", width=8, height=8)
ggplot(DE_homolog2, aes(y=kin20downdiffa2,x=genomiccomp, fill=DE_status)) + 
  scale_fill_manual(values = c("firebrick2","light grey","dodgerblue2"),labels=c("A2-biased","Not-biased","A1-biased"), name="Expression") +
  geom_boxplot(outlier.shape=NA) +
  ylim(-6,6) +  
  labs(x='Genomic compartment') + 
  labs(y='Difference of TE number [10-20k downstream] (A1-A2)') +
  theme(axis.title.x = element_text(size=10,colour = "black"),axis.title.y = element_text(size=10,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=10),axis.text.y = element_text(colour="black",size=10))
dev.off()

