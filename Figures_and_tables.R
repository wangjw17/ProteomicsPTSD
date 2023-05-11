##############################################################
### FIGURE MAKING IN PROTEOMICS/MIRNA PROJECT
### Jiawei Wang, 04/12/2022 
##############################################################
library(dplyr)
library(ggplot2)

#### miRNA Clustering dendrogram with trait ####
rm(list=ls())
load('rawdata/dat_upmc_mirna_all1.5.RData')
getFKPM <- function(dat, log2=TRUE, minute=1e-3, trim=5, scale=1e6, suffix=TRUE, datLen){
  sel <- (trim+1):dim(dat)[2]
  seldat <- dat[,sel]
  sums <- colSums(seldat)/scale
  seldat <- t(t(seldat)/sums)
  len <- datLen/1000
  seldat <- seldat/len
  if(log2)  seldat <- log(seldat+minute, base=2)
  if(suffix) seldat <- cbind(dat[,1:trim], seldat)
  return(seldat)
}

## all transcripts
datFPKM <- getFKPM(datExpr, trim=0, minute=1, suffix=F, datLen = datProbes$Length)
sel.mir <- which(datProbes$TranscriptType=="miRNA")
sel.sam <- which(datMeta$PrimaryDx!="MDD")
datFPKM <- datFPKM[sel.mir,sel.sam]
datProbes1 <- datProbes[sel.mir,]
datExpr1 <- datExpr[sel.mir,sel.sam]
datMeta1 <- datMeta[sel.sam,]

library(WGCNA)
sampleTree2 = hclust(dist(t(datExpr1)), method = "average")
datTraits <- datMeta1[c('PrimaryDx','region','Sex','AgeDeath','Race','RIN','PMI')]
names(datTraits)[c(2,5)] <- c("Brain region","Ancestry")
datTraits$PrimaryDx <- as.numeric(datTraits$PrimaryDx)
datTraits$`Brain region` <- as.numeric(datTraits$`Brain region`)
datTraits$Ancestry <- as.numeric(datTraits$Ancestry)
datTraits$Sex <- as.numeric(datTraits$Sex)
traitColors = numbers2colors(datTraits, signed = F, colors = );
# Plot the sample dendrogram and the colors underneath.
pdf('manuscript/figures_v5/1c_dendro.pdf', width=12, height=4)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "miRNA sample dendrogram and trait heatmap")
dev.off()


### Heatmap of top miRNAs 
sel.mir2 <- order(apply(datFPKM,1,var), decreasing = T)[1:50]
datExpr2 <- datFPKM[sel.mir2,] %>% t() %>% scale %>% t() %>% as.data.frame
datProbes2 <- datProbes1[sel.mir2,]
df <- reshape2::melt(cbind(datProbes2[,1],datExpr2))
df$Sample <- factor(df$Sample, levels=names(datExpr2)[sampleTree2$order])
names(df) <- c("miRNA",'Sample','Expr')
pdf('manuscript/figures_v5/1d_heatmap_unf.pdf', width=12, height=4)
ggplot(df, aes(Sample, miRNA, fill= Expr)) + 
  geom_tile(colour="white", size=.2) +
  scale_fill_gradient2(low = "blue4", high = "red3", mid = "white",# midpoint = 0, limit = c(-4,4), 
                       space = "Lab", name="Relative\nAbundance") +
  # scale_y_discrete(position = "right") +
  labs(x=NULL,y=NULL)+
  theme(#axis.text.y.right = element_text(face="italic"),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    # axis.ticks = element_blank(), 
    panel.background=element_blank(), legend.position = "none")
dev.off()

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "RdBu"))(25) ##Accent, Set2, Set3, RdYlBu(-),Spectral(-)
heatmap(as.matrix(datExpr2[-1]), #ColSideColors = datMeta$color[datMeta$Region=="dlPFC"],
        # Colv = NA, 
        scale="none", col=coul, labCol = NA, cexRow = 0.6)


#### miRNA Volvano plot in PTSD ####
rm(list=ls())
## dlPFC+sgPFC in ggplot
library(ggplot2)
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
sig01$Dx <- "MDD"; sig02$Dx <- "PTSD"; sig.dl <- rbind(sig01, sig02); sig.dl$BA <- 'dlPFC'
load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
sig01$Dx <- "MDD"; sig02$Dx <- "PTSD"; sig.sg <- rbind(sig01, sig02); sig.sg$BA <- 'sgPFC'
sig <- rbind(sig.dl, sig.sg)
sig$dir <- 0
sig$dir[sig$log2FoldChange>log(1.5)/log(2) & sig$pvalue<.05] <- 1
sig$dir[sig$log2FoldChange< -log(1.5)/log(2) & sig$pvalue<.05] <- -1
sig$dir <- factor(sig$dir, levels=c(0,1,-1))
sig$Transcriptid <- as.character(sig$Transcriptid)
sig$label <- NA; sig$label[sig$dir!=0] <- sig$Transcriptid[sig$dir!=0]
sig$label <- sig$label %>% tolower %>% gsub('mir','hsa-miR-',.) %>% gsub('miR-let','Let-',.)
sig$Dx <- factor(sig$Dx, levels=c("PTSD","MDD"))
sig$dir1 <- 0
sig$dir1[sig$log2FoldChange > 0 & sig$pvalue < .05] <- 1
sig$dir1[sig$log2FoldChange < 0 & sig$pvalue < .05] <- -1
sig$dir1[sig$dir1==1 & sig$dir==0] <- 0.5
sig$dir1[sig$dir1==-1 & sig$dir==0] <- -0.5
sig$dir1 <- factor(sig$dir1, levels=c(0,1,-1,.5,-.5))
pdf('manuscript/figures_v5/1a_volcano.pdf', width=7, height=8)
ggplot(subset(sig,Dx=="PTSD"), aes(x=log2FoldChange, y=-log10(pvalue), label=label)) +
  geom_point(aes(fill=dir1), shape = 21, colour = "black", size = 4) +
  # scale_fill_manual(values = c('grey80','firebrick','steelblue','coral','skyblue')) +
  scale_fill_manual(values = c('grey80','firebrick','forestgreen','coral','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*log(1.5)/log(2), lty=2, col="blue") +
  # geom_text(fontface = "bold",position=position_jitter(width=.1,height=.1)) +
  ggrepel::geom_text_repel()+#fontface = "bold", size=4) +
  # facet_grid(rows = vars(Dx), cols = vars(BA), switch="y") + 
  facet_grid(BA~., switch="y") +
  xlim(-3,3) + xlab("log2(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")
dev.off()

### simplified for schematic purpose
ggplot(subset(sig,Dx=="PTSD"), aes(x=log2FoldChange, y=-log10(pvalue), label=label)) +
  geom_point(aes(fill=dir1), shape = 21, colour = "black", size = 4) +
  # scale_fill_manual(values = c('grey80','firebrick','steelblue','coral','skyblue')) +
  scale_fill_manual(values = c('grey80','firebrick','forestgreen','coral','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*log(1.5)/log(2), lty=2, col="blue") +
  # geom_text(fontface = "bold",position=position_jitter(width=.1,height=.1)) +
  # ggrepel::geom_text_repel()+#fontface = "bold", size=4) +
  # facet_grid(rows = vars(Dx), cols = vars(BA), switch="y") + 
  # facet_grid(BA~., switch="y") +
  xlim(-3,3) + xlab("log2(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")


#### miRNA WGCNA dendrogram ####
rm(list=ls())
load('rawdata/dat_upmc_mirna_all1.5.RData')
getFKPM <- function(dat, log2=TRUE, minute=1e-3, trim=5, scale=1e6, suffix=TRUE, datLen){
  sel <- (trim+1):dim(dat)[2]
  seldat <- dat[,sel]
  sums <- colSums(seldat)/scale
  seldat <- t(t(seldat)/sums)
  len <- datLen/1000
  seldat <- seldat/len
  if(log2)  seldat <- log(seldat+minute, base=2)
  if(suffix) seldat <- cbind(dat[,1:trim], seldat)
  return(seldat)
}

## all transcripts
datFPKM <- getFKPM(datExpr, trim=0, minute=1, suffix=F, datLen = datProbes$Length)
sel.mir <- which(datProbes$TranscriptType=="miRNA")
datFPKM <- datFPKM[sel.mir,]
datProbes1 <- datProbes[sel.mir,]
datExpr1 <- datExpr[sel.mir,]

sel0 <- !is.na(datProbes1$Genename) & rowMeans(datFPKM)>0.5
# sam0 <- datMeta$region == "sgPFC" ## dlPFC
sam0 <- 1:nrow(datMeta)
# sam0 <- datMeta$PrimaryDx!="MDD"
# sam0 <- datMeta$PrimaryDx!="PTSD"
datFPKM <- datFPKM[sel0,sam0]
datMeta <- datMeta[sam0,]
datProbes2 <- datProbes1[sel0,]
datExpr2 <- datExpr1[sel0,sam0]

### look at batch ###
batch1 <- paste0("RTMIR0",c(40:62)) ## 12
# batch2 <- "RTMIR098 RTMIR099 RTMIR093 RTMIR090 RTMIR096 RTMIR097 RTMIR095 RTMIR085 RTMIR084 RTMIR087 RTMIR081 RTMIR083 RTMIR088 RTMIR066 RTMIR080 RTMIR092 RTMIR075 RTMIR077 RTMIR070 RTMIR071 RTMIR072 RTMIR079 RTMIR065 RTMIR063 RTMIR069 RTMIR101 RTMIR082 RTMIR086 RTMIR078 RTMIR104 RTMIR100 RTMIR103" %>% strsplit(split=" ") %>% unlist
batch2 <- paste0("RTMIR",c(paste0("0",63:99),100:104)) ## 36
batch3 <- paste0("RTMIR", c(120:132)) ## 12
batch4 <- paste0("RTMIR0",c(paste0("0",1:9),10:39)) ## 36
batch5 <- paste0("RTMIR",105:119) ## 12
datMeta$batch[datMeta$Sample.mir %in% batch1] <- 1
datMeta$batch[datMeta$Sample.mir %in% batch2] <- 2
datMeta$batch[datMeta$Sample.mir %in% batch3] <- 3
datMeta$batch[datMeta$Sample.mir %in% batch4] <- 4
datMeta$batch[datMeta$Sample.mir %in% batch5] <- 5

### combat correction
library(sva)
# load('DIA_new.RData')
mod = model.matrix(~PrimaryDx+AgeDeath+Sex+RIN+PMI, data=datMeta)
# datFPKM.combat = ComBat(dat = as.matrix(datFPKM), batch = factor(datMeta$region), mod = mod)
datFPKM.combat = ComBat(dat = as.matrix(datFPKM), batch = factor(datMeta$batch), mod = mod)
datFPKM.preCombat = datFPKM
datFPKM = as.data.frame(datFPKM.combat)

### sample selection after batch correction
sam1 <- datMeta$PrimaryDx!="MDD"
datFPKM <- datFPKM[,sam1]
datMeta <- datMeta[sam1,]
datExpr2 <- datExpr1[,sam1]

# ## remove covariates
dat <- datFPKM %>% t()
# datMeta0 <- datMeta[c('AgeDeath', 'RIN', 'PMI', 'Sex')]
# datMeta0$Sex <- as.numeric(datMeta0$Sex)
# lm.cov <- lm(dat ~ AgeDeath + RIN + PMI + Sex, data=datMeta0) ##not including Race
# dat.ori <- dat
# dat <- dat - as.matrix(datMeta0[c('AgeDeath', 'RIN', 'PMI', 'Sex')]) %*% lm.cov$coefficients[2:5,]
power <- 6
library(WGCNA)
net <- blockwiseModules(datExpr = dat, power = power,
                        TOMType = "signed", minModuleSize = 10,
                        reassignThreshold = 0.2, mergeCutHeight = 0,
                        numericLabels = F, pamRespectsDendro = F,
                        saveTOMs = F, 
                        saveTOMFileBase = "catTOM",
                        verbose = 3)
# save(net, datExpr2, file="result_mirna/WGCNA_PTSD_dl_sg.RData")
pdf('manuscript/figures_v5/2a_wgcna_dendrogram_traits.pdf', width=8, height=5)
datTraits <- datMeta[c('PrimaryDx','region','Sex','AgeDeath','Race','RIN','PMI')]
names(datTraits)[c(2,5)] <- c("Brain region","Ancestry")
datTraits$PrimaryDx <- as.numeric(datTraits$PrimaryDx)
datTraits$`Brain region` <- as.numeric(datTraits$`Brain region`)
# datTraits$Ancestry <- as.numeric(datTraits$Ancestry)
datTraits$Race.AA <- as.numeric(datTraits$Ancestry=="AA")
datTraits$Race.EA <- as.numeric(datTraits$Ancestry=="CAUC")
datTraits$Sex <- as.numeric(datTraits$Sex)
datTraits <- datTraits[,-5]
datTraits.cor <- cor(datTraits, t(datFPKM)) %>% t()
traitColors = numbers2colors(datTraits.cor, signed = T)
traitColors <- cbind(net$colors[net$blockGenes[[1]]], traitColors)
plotDendroAndColors(net$dendrograms[[1]], traitColors,
                    main = "miRNA dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    groupLabels = c('Module', names(datTraits)), 
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#### miRNA module-trait heatmap ####
library(magrittr)
library(stats)
library(WGCNA)
library(flashClust)
library(biomaRt)
library(corrplot)
library(ggplot2)
# library(cqn)
# library(sva)
# library(limma)
library(statmod)
library(dplyr)
# library(pSI)
# library(pSI.data)
library(gplots)
library(stringr)
# options(stringsAsFactors = FALSE, digits = 3)
# theme_update(plot.title = element_text(hjust = 0.5))
rm(list=ls())
load('result_2022/WGCNA_CP_ds_combatBatch.RData')

## dendrogram
# pdf('result_0703//wgcna_dendrogram_amp_ca.pdf', width=12)
# plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
#                     main = "Single block gene dendrogram and module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)


## module-trait cor
datMeta1 <- datMeta[c('Sex','AgeDeath','RIN','PMI')]
# datMeta1$MDD <- as.numeric(datMeta$PrimaryDx=="MDD")
datMeta1$PTSD <- as.numeric(datMeta$PrimaryDx=="PTSD")
datMeta1$Sex <- as.numeric(datMeta1$Sex)
datMeta1$Race.AA <- as.numeric(datMeta$Race=="AA")
datMeta1$Race.EA <- as.numeric(datMeta$Race=="CAUC")
# moduleColors <- labels2colors(net$colors)
nGenes = ncol(dat);
# MEs0 = moduleEigengenes(dat, net$colors)$eigengenes
MEs = orderMEs(MEs0)
names(MEs) <- names(MEs) %>% gsub('ME','',.) %>% as.numeric %>% labels2colors
nSamples <- dim(datMeta1)[1]
moduleTraitCor = cor(MEs, datMeta1, use = "p");# Pearson correlation
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf('manuscript/figures_v5/2_wgcna_traitcor_ds_cp.pdf', width=8, height=4)
par(mar = c(5, 5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datMeta1),
               yLabels = names(MEs),
               ySymbols = names(MEs) %>% gsub("ME","",.),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()




#### miRNA module membership vs GS ####
load('result_2022/WGCNA_CP_ds_combatBatch.RData')
plot(data=subset(all,module==1), trait.GS.PTSD~DxmoduleMembership)
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
all1 <- subset(all, module==1)
names(all1)[c(9,12)] <- c("y","x")
library(ggpmisc)
pdf('manuscript/figures_v5/2c_regression.pdf', height=5, width=8)
ggplot(all1, aes(x=x,y=y)) +
  geom_point(pch=16, size=3) + xlab('Module Membership') + ylab('PTSD correlation')+
  geom_smooth(method='lm',formula=y~x) +
  # geom_text(x=-.95, y=-0.15, label = lm_eqn(all1), parse = TRUE)+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_bw()
dev.off()


#### miRNA ME1 volcano plot ####
rm(list=ls())
library(ggplot2)
library(ggrepel)
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
load('result_2022/WGCNA_CP_ds_combatBatch.RData')
mod <- "turquoise"
alpha_fc = log(1.5)/log(2)
# de_mod <- de_sp[de_sp$mod.CP==mod,]
# de_mod$GN.sig <- de_mod$GN; de_mod$GN.sig[de_mod$PTSD.P.Value>.05] <- NA
sig02$mod <- all$module[match(sig02$Geneid,all$Protein)]
de_mod <- subset(sig02, mod==1)
de_mod$GN.sig <- de_mod$Transcriptid; de_mod$GN.sig[de_mod$pvalue>.05] <- NA
de_mod$GN.sig <- de_mod$GN.sig %>% tolower %>% gsub('mir','hsa-miR-',.)
de_mod$dir <- 0
de_mod$dir[de_mod$log2FoldChange>0 & de_mod$pvalue<.05] <- 1
de_mod$dir[de_mod$log2FoldChange<0 & de_mod$pvalue<.05] <- -1
de_mod$dir[de_mod$log2FoldChange>0 & de_mod$log2FoldChange<alpha_fc & de_mod$pvalue<.05] <- .5
de_mod$dir[de_mod$log2FoldChange<0 & de_mod$log2FoldChange>-alpha_fc & de_mod$pvalue<.05] <- -.5
de_mod$dir <- factor(de_mod$dir, levels=c(0,1,-1,.5,-.5))
pdf('manuscript/figures_v5/2d_volcano1.pdf', height=5)
# ggplot(de_mod, aes(x=log2FoldChange,y=-log10(pvalue),fill=dir)) +
#   # geom_point(alpha=.3+0.7*(de_mod$pvalue<.05)) +
#   geom_point(shape = 21, colour = "black", size = 4) +
#   scale_fill_manual(values = c('grey80','firebrick','forestgreen')) +
#   geom_hline(yintercept=log10(20), linetype="dashed", color = "red") +
#   geom_text_repel(aes(label=GN.sig), size = 5,
#                   box.padding = unit(0.35, "lines"),
#                   point.padding = unit(0.3, "lines")) +
#   xlab('log10FC') + ylab('-log10(Pvalue)') + labs(title='Module turquoise') +
#   theme_bw()
ggplot(de_mod, aes(x=log2FoldChange, y=-log10(pvalue), label=GN.sig)) +
  geom_point(aes(fill=dir), shape = 21, colour = "black", size = 4) +
  scale_fill_manual(values = c('grey80','firebrick','forestgreen','coral','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*log(1.5)/log(2), lty=2, col="blue") +
  ggrepel::geom_text_repel()+#fontface = "bold", size=4) +
  # facet_grid(BA~., switch="y") +
  # xlim(-3,3) + 
  xlab("log2(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")
dev.off()

#### protein PTSD DEP volcano & venn ####
library(ggVennDiagram)
rm(list=ls())
# dep.dl <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
# dep.sg <- read.csv('result_0518/DE_MDD_PTSD_sg.csv')
dep.dl <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
dep.sg <- read.csv('result_0518/DE_MDD_PTSD_sg_incA_new1.csv')
dep <- merge(dep.dl, dep.sg, by="EN")
names(dep) <- names(dep) %>% gsub('x','dl',.) %>% gsub('y','sg',.)
# x <- list(`PTSD dlPFC`=dep$EN[dep$PTSD.P.Value.dl<.05],
#           `PTSD sgPFC`=dep$EN[dep$PTSD.P.Value.sg<.05],
#           `MDD dlPFC`=dep$EN[dep$MDD.P.Value.dl<.05],
#           `MDD sgPFC`=dep$EN[dep$MDD.P.Value.sg<.05])
x <- list(`PTSD dlPFC`=dep$EN[dep$PTSD.P.Value.dl<.05],
          `PTSD sgPFC`=dep$EN[dep$PTSD.P.Value.sg<.05])
pdf('manuscript/figures_v5/3b_vennProtein.pdf', width=6, height=3)
ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey80", high = "coral3")
dev.off()

### style as in NN21
c3 <- dep[c("PTSD.P.Value.dl","PTSD.P.Value.sg")]
names(c3) <- c("PTSD dlPFC","PTSD sgPFC")
c3 <- (c3<.05)
library(limma)
a <- vennCounts(c3)
pdf('manuscript/figures_v5/3b_vennProtein.pdf')
vennDiagram(a, circle.col = 2:3)
dev.off()

library(ggplot2)
rm(list=ls())
# sig01 <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
# sig02 <- read.csv('result_0518/DE_MDD_PTSD_sg.csv')
sig01 <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
sig02 <- read.csv('result_0518/DE_MDD_PTSD_sg_incA_new1.csv')
sig01$BA <- "dlPFC"; sig02$BA <- "sgPFC"
names(sig01) <- names(sig01) %>% gsub("MDD.|PTSD.","",.)
names(sig02) <- names(sig02) %>% gsub("MDD.|PTSD.","",.)
sig <- rbind(sig01[c(2:9,18:21)], sig01[10:21], sig02[c(2:9,18:21)], sig02[10:21])
sig$dir <- 0
thres.logfc <- log(1.5)/log(10)
sig$dir[sig$logFC>thres.logfc & sig$P.Value<.05] <- 1
sig$dir[sig$logFC< -thres.logfc & sig$P.Value<.05] <- -1
sig$dir <- factor(sig$dir, levels=c(0,1,-1))
sig$GN <- as.character(sig$GN)
sig$label <- NA; sig$label[sig$dir!=0] <- sig$GN[sig$dir!=0]
sig$Dx <- c(rep(c("MDD","PTSD"), each=nrow(sig01)), rep(c("MDD","PTSD"), each=nrow(sig02))) %>%  factor(levels=c("PTSD","MDD"))
sig$dir1 <- 0
sig$dir1[sig$logFC > 0 & sig$P.Value < .05] <- 1
sig$dir1[sig$logFC < 0 & sig$P.Value < .05] <- -1
sig$dir1[sig$dir1==1 & sig$dir==0] <- 0.5
sig$dir1[sig$dir1==-1 & sig$dir==0] <- -0.5
sig$dir1 <- factor(sig$dir1, levels=c(0,1,-1,.5,-.5))
pdf('manuscript/figures_v5/3a_volcano.pdf', width=5, height=6)
ggplot(subset(sig,Dx=="PTSD"), aes(x=logFC, y=-log10(P.Value), label=label)) +
  geom_point(aes(fill=dir1), shape = 21, colour = "black", size = 4) +
  # scale_color_manual(values = c('grey','firebrick','steelblue','coral','turquoise')) +
  scale_fill_manual(values = c('grey80','firebrick','forestgreen','coral','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*thres.logfc, lty=2, col="blue") +
  # geom_text(fontface = "bold",position=position_jitter(width=.1,height=.1)) +
  ggrepel::geom_text_repel(fontface = "bold", size=4) +
  # facet_grid(rows=vars(Dx), cols=vars(BA), switch="y") + 
  facet_grid(BA~., switch="y") +
  xlim(-.4,.4) + xlab("log10(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")
dev.off()




#### protein coverage of transcript types ####
#### miRNA enrichment of DEPs, PTSD ####
rm(list=ls())
df.p <- read.csv('result_mirna/pair_dlpfc_ptsd_share_511_enrichment_new1.csv'); df.p$Dx <- "PTSD"
df.p$miRNA <- df.p$miRNA %>% tolower %>% gsub('mir','hsa-miR-',.) ## NOMENCLATURE
df.p$miRNA <- factor(df.p$miRNA, levels=df.p$miRNA[order(df.p$p.chisq)] %>% rev)
df.p <- subset(df.p, miRNA!="")
df.p <- subset(df.p, miRNA %in% df.p$miRNA[df.p$p.chisq<.05])
pdf('manuscript/figures_v5/3d1_enrichment_dl1.pdf', width=5, height=4)
ggplot(df.p, aes(miRNA, -log10(p.chisq)*p.sign)) +
  geom_bar(stat="identity", position="dodge", color="black", fill="grey80", width = .6) +
  geom_hline(yintercept=c(-log10(0.05)), color="grey60", linetype="dashed") +
  scale_fill_gradient2(low="blue", mid="white", high="firebrick2") +
  xlab('miRNAs') + ylab('signed -log10(P)') + ylim(0, 10) +
  theme_bw() +
  # theme_classic() +
  labs(title = "miRNA enrichment of DEP pairs, dlPFC") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8)) + coord_flip()
dev.off()
# df.p$BA <- "dlPFC"
# df.p1 <- df.p

pdf('manuscript/figures_v5/3d2_enrichment_sg1.pdf', width=5, height=1.6)
df.p <- read.csv('result_mirna/pair_sgpfc_ptsd_share_511_enrichment_new1.csv'); df.p$Dx <- "PTSD"
df.p$miRNA <- df.p$miRNA %>% tolower %>% gsub('mir','hsa-miR-',.) ## NOMENCLATURE
df.p$miRNA <- factor(df.p$miRNA, levels=df.p$miRNA[order(df.p$p.chisq)] %>% rev)
df.p <- subset(df.p, miRNA!="")
df.p <- subset(df.p, miRNA %in% df.p$miRNA[df.p$p.chisq<.05])
df.p$BA <- "sgPFC"
# df.p <- rbind(df.p1, df.p)
ggplot(df.p, aes(miRNA, -log10(p.chisq)*p.sign)) + 
  geom_bar(stat="identity", position="dodge", color="black", fill="grey80", width = .6) +
  geom_hline(yintercept=c(-log10(0.05)), color="grey60", linetype="dashed") +
  # scale_fill_gradient2(low="blue", mid="white", high="firebrick2") +
  xlab('miRNAs') + ylab('signed -log10(P)') + ylim(0, 10) +
  theme_bw() +
  # theme_classic() +
  # facet_grid(BA~., switch="y") + 
  labs(title = "miRNA enrichment of DEP pairs, sgPFC") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8)) + coord_flip()
dev.off()



#### logFC comparison between mRNA & protein ####
rm(list=ls())
# load('result_mirna/pair_dlpfc_CM.RData')
load('result_mirna/pair_dlpfc_CP.RData')
load('../transcriptomics_ctx/dat_ulval/deseq_A9_m3_ulval_keep3.RData')
# load('dat_upmc/deseq_A9_m3_upmc_keep3.RData')
# load('dat_lval/deseq_A9_m3_lval_keep3.RData')
dep <- read.csv('result_200518/DE_MDD_PTSD_dl_incA_new1.csv'); names(dep)[19] <- "Genename"
dep <- subset(dep, Genename %in% df$Var2)
degs <- merge(dep, sig01, by="Genename")
degs <- merge(degs, sig02, by="Genename")
names(degs) <- names(degs) %>% gsub(".x",".MDD",.) %>% gsub('.y','.PTSD',.)

pdf('manuscript/figures/comparison_mRNA-protein_mdd.pdf', width=6, height=4.5)
sel <- apply(degs[c('PTSD.P.Value','pvalue.PTSD')], 1, max) < .05 ## both moderately significant
degs$`P-value` <- apply(degs[c('PTSD.P.Value','pvalue.PTSD')],1,min)
degs$miRNA="No miRNA associated"; #degs$miRNA[degs$Genename %in% c("CD59", "DPP6")] <- "miRNA associated"
ggplot(degs[sel,], aes(x=PTSD.logFC, y=log2FoldChange.PTSD, label=Genename))+#, color=factor(miRNA))) +
  geom_point(aes(size=-log10(`P-value`)), alpha=0.8) +
  ggrepel::geom_text_repel(size=4) +
  scale_color_manual(values = c("red","black")) +
  geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
  xlab('Protein log10FC') + ylab('RNA log2FC') + labs(title="PTSD dlPFC") +
  theme_bw()
sel <- apply(degs[c('MDD.P.Value','pvalue.MDD')], 1, max) < .05 ## both moderately significant
degs$`P-value` <- apply(degs[c('MDD.P.Value','pvalue.MDD')],1,min)
degs$miRNA="No miRNA associated";# degs$miRNA[degs$Genename %in% c("CD59", "DPP6")] <- "miRNA associated"
ggplot(degs[sel,], aes(x=MDD.logFC, y=log2FoldChange.MDD, label=Genename))+#, color=factor(miRNA))) +
  geom_point(aes(size=-log10(`P-value`)), alpha=0.8) +
  ggrepel::geom_text_repel(size=4) +
  scale_color_manual(values = c("red","black")) +
  geom_hline(yintercept=0, lty=2) + geom_vline(xintercept=0, lty=2) +
  xlab('Protein log10FC') + ylab('RNA log2FC') + labs(title="MDD dlPFC") +
  theme_bw()
dev.off()

#### protein module-trait heatmap ####
# pdf('manuscript/figures/3_moduleTrait.pdf', width=8, height=8)
rm(list=ls())
library(WGCNA)
library(dplyr)
## PTSD, dlPFC
load('result_0703/WGCNA_CP_dl.RData')
nSamples <- dim(datMeta)[1]
datMeta$PTSD <- as.numeric(datMeta$PrimaryDx=="PTSD")
datMeta$Ancestry <- datMeta$Race %>% gsub("Asian",0,.) %>% gsub('B',1,.) %>% gsub('W',-1,.) %>% as.numeric
datMeta1 <- datMeta[c("PTSD","AgeDeath","PMI","Sex","Ancestry")]
datMeta1$Sex <- as.numeric(datMeta1$Sex)
MEs <- MEs[,order(moduleTraitCor[,'ptsd'], decreasing = T)]
moduleTraitCor = cor(MEs, datMeta1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# textMatrix = paste(signif(moduleTraitCor, 2), " (",
#                    signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = signif(-log10(moduleTraitPvalue), 2)
dim(textMatrix) = dim(moduleTraitCor)
id.grey <- which(names(MEs)=="grey")
pdf('manuscript/figures_v5/4a_module_trait_CP_dl.pdf')
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor[-id.grey,],
               xLabels = names(datMeta1),
               yLabels = names(MEs)[-id.grey],
               ySymbols = names(MEs)[-id.grey],
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix[-id.grey,],
               setStdMargins = FALSE,
               cex.text = 0.8,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
## MDD, dlPFC
rm(list=ls())
load('result_0703/WGCNA_CM_dl.RData')
names(MEs) <- names(MEs) %>% gsub('ME','',.)
nSamples <- dim(datMeta)[1]
datMeta$MDD <- as.numeric(datMeta$PrimaryDx=="MDD")
datMeta$Ancestry <- datMeta$Race %>% gsub("Asian",0,.) %>% gsub('B',1,.) %>% gsub('W',-1,.) %>% as.numeric
datMeta1 <- datMeta[c("MDD","AgeDeath","PMI","Sex","Ancestry")]
datMeta1$Sex <- as.numeric(datMeta1$Sex)
MEs <- MEs[,order(moduleTraitCor[,'mdd'], decreasing = T)]
moduleTraitCor = cor(MEs, datMeta1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# textMatrix = paste(signif(moduleTraitCor, 2), " (",
#                    signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = signif(-log10(moduleTraitPvalue), 2)
dim(textMatrix) = dim(moduleTraitCor)
id.grey <- which(names(MEs)=="grey")
# pdf('result_0128/module_trait_CM_dl.pdf')
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor[-id.grey,],
               xLabels = names(datMeta1),
               yLabels = names(MEs)[-id.grey],
               ySymbols = names(MEs)[-id.grey],
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix[-id.grey,],
               setStdMargins = FALSE,
               cex.text = 0.8,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#### protein modular volcano plots: CP-skyblue/red, CM-darkred/grey60/black ####
rm(list=ls())
library(ggplot2)
library(ggrepel)
# dep <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
# load('result_0703/WGCNA_CP_dl.RData')
# dep$mod.CP <- all$module[match(dep$EN, all$Protein)]
# load('result_0703/WGCNA_CM_dl.RData')
# dep$mod.CM <- all$module[match(dep$EN, all$Protein)]
# write.csv(dep[-1], file="result_0128/DE_MDD_PTSD_dl_mod_incA_new1.csv", row.names=F)

pdf('manuscript/figures_v5/4b_volcanoMod2.pdf', height=4, width=4)
load('DIA_new.RData')
# de_sp <- read.csv('result_0128/DE_MDD_PTSD_dl_mod.csv')
de_sp <- read.csv('result_0128/DE_MDD_PTSD_dl_mod_incA_new1.csv')
# de_sp <- de_sp[de_sp$EN %in% datProbes$EN[datProbes$DB=="sp"],]

mod <- "skyblue" ## for PTSD, red/skyblue/tan/royalblue
de_mod <- de_sp[de_sp$mod.CP==mod,]
de_mod$GN.sig <- de_mod$GN; de_mod$GN.sig[de_mod$PTSD.P.Value>.05] <- NA
thres.logfc = log10(1.5)
de_mod$dir <- 1
de_mod$dir[de_mod$PTSD.logFC>0 & de_mod$PTSD.P.Value<.05] <- 4
de_mod$dir[de_mod$PTSD.logFC>thres.logfc & de_mod$PTSD.P.Value<.05] <- 2
de_mod$dir[de_mod$PTSD.logFC<0 & de_mod$PTSD.P.Value<.05] <- 5
de_mod$dir[de_mod$PTSD.logFC< -thres.logfc & de_mod$PTSD.P.Value<.05] <- 3
de_mod$dir <- factor(de_mod$dir, levels=1:5)
ggplot(de_mod, aes(x=PTSD.logFC, y=-log10(PTSD.P.Value))) +
  geom_point(aes(fill=dir), shape = 21, colour = "black", size = 4) +
  # scale_color_manual(values = c('grey','firebrick','steelblue','coral','turquoise')) +
  scale_fill_manual(values = c('grey80','coral','lightgreen')) +
  # scale_fill_manual(values = c('grey80','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*thres.logfc, lty=2, col="blue") +
  geom_text_repel(aes(label=GN.sig), size = 5,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")) +
  xlim(-.25,.25) + xlab("log10(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")

mod <- "red" ## for PTSD, red/skyblue/tan/royalblue
de_mod <- de_sp[de_sp$mod.CP==mod,]
de_mod$GN.sig <- de_mod$GN; de_mod$GN.sig[de_mod$PTSD.P.Value>.05] <- NA
thres.logfc = log10(1.5)
de_mod$dir <- 1
de_mod$dir[de_mod$PTSD.logFC>0 & de_mod$PTSD.P.Value<.05] <- 4
de_mod$dir[de_mod$PTSD.logFC>thres.logfc & de_mod$PTSD.P.Value<.05] <- 2
de_mod$dir[de_mod$PTSD.logFC<0 & de_mod$PTSD.P.Value<.05] <- 5
de_mod$dir[de_mod$PTSD.logFC< -thres.logfc & de_mod$PTSD.P.Value<.05] <- 3
de_mod$dir <- factor(de_mod$dir, levels=1:5)
ggplot(de_mod, aes(x=PTSD.logFC, y=-log10(PTSD.P.Value))) +
  geom_point(aes(fill=dir), shape = 21, colour = "black", size = 4) +
  # scale_color_manual(values = c('grey','firebrick','steelblue','coral','turquoise')) +
  scale_fill_manual(values = c('grey80','coral','lightgreen')) +
  # scale_fill_manual(values = c('grey80','lightgreen')) +
  geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
  geom_vline(xintercept = c(-1,1)*thres.logfc, lty=2, col="blue") +
  geom_text_repel(aes(label=GN.sig), size = 5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  xlim(-.25,.25) + xlab("log10(Fold Change)") + ylab("-log10(P value)") + 
  theme_bw() + theme(legend.position = "none")

# mod <- "grey60" ## for MDD, darkred/grey60/turquoise
# de_mod <- de_sp[de_sp$mod.CM==mod,]
# de_mod$GN.sig <- de_mod$GN; de_mod$GN.sig[de_mod$MDD.P.Value>.05] <- NA
# thres.logfc = log10(1.5)
# de_mod$dir <- 1
# de_mod$dir[de_mod$MDD.logFC>0 & de_mod$MDD.P.Value<.05] <- 4
# de_mod$dir[de_mod$MDD.logFC>thres.logfc & de_mod$MDD.P.Value<.05] <- 2
# de_mod$dir[de_mod$MDD.logFC<0 & de_mod$MDD.P.Value<.05] <- 5
# de_mod$dir[de_mod$MDD.logFC< -thres.logfc & de_mod$MDD.P.Value<.05] <- 3
# de_mod$dir <- factor(de_mod$dir, levels=1:5)
# 
# ggplot(de_mod, aes(x=MDD.logFC, y=-log10(MDD.P.Value))) +
#   geom_point(aes(fill=dir), shape = 21, colour = "black", size = 4) +
#   # scale_color_manual(values = c('grey','firebrick','steelblue','coral','turquoise')) +
#   scale_fill_manual(1:5,values = c('grey80','firebrick','lightgreen','coral')) +
#   # scale_fill_manual(values = c('grey80','lightgreen')) +
#   geom_hline(yintercept = -log10(0.05), lty=2, col="red") +
#   geom_vline(xintercept = c(-1,1)*thres.logfc, lty=2, col="blue") +
#   geom_text_repel(aes(label=GN.sig), size = 5,
#                   box.padding = unit(0.35, "lines"),
#                   point.padding = unit(0.3, "lines")) +
#   xlim(-.25,.25) + xlab("log10(Fold Change)") + ylab("-log10(P value)") +
#   theme_bw() + theme(legend.position = "none")
dev.off()

#### protein modular eigengene plots: CP-skyblue/red, CM-darkred/grey60/black ####
rm(list=ls())
# pdf('result_0128/module_eigengene_dl.pdf')
# pdf('manuscript/figures/3_eigenprot.pdf', width=4, height=4)
pdf('manuscript/figures_v5/4b_eigenprot2.pdf', width=4, height=4)
# load('result_0703/WGCNA_CM_dl.RData')
# df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$MEgrey60)
# ggplot(df, aes(x=PrimaryDx, y=Eigengene)) +
#   geom_boxplot() +
#   # scale_fill_brewer(palette="Blues") +
#   # scale_fill_manual(values=c("grey60","grey80")) +
#   geom_jitter(aes(color=PrimaryDx), size=4, shape=16, position=position_jitter(0.1)) +
#   labs(title="MDD module grey60") + theme_classic()
# df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$MEdarkred)
# ggplot(df, aes(x=PrimaryDx, y=Eigengene)) +
#   geom_boxplot() +
#   # scale_fill_brewer(palette="Blues") +
#   # scale_fill_manual(values=c("grey60","grey80")) +
#   geom_jitter(aes(color=PrimaryDx), size=4, shape=16, position=position_jitter(0.1)) +
#   labs(title="MDD module darkred") + theme_classic()
# df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$MEblack)
# ggplot(df, aes(x=PrimaryDx, y=Eigengene)) +
#   geom_boxplot() +
#   # scale_fill_brewer(palette="Blues") +
#   # scale_fill_manual(values=c("grey60","grey80")) +
#   geom_jitter(aes(color=PrimaryDx), size=4, shape=16, position=position_jitter(0.1)) +
#   labs(title="MDD module black") + theme_classic()
load('result_0703/WGCNA_CP_dl.RData')
df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$red)
ggplot(df, aes(x=PrimaryDx, y=Eigengene)) +
  geom_boxplot() +
  # scale_fill_brewer(palette="Blues") +
  # scale_fill_manual(values=c("grey60","grey80")) +
  geom_jitter(aes(fill=PrimaryDx), color="black", size=4, shape=21, position=position_jitter(0.1)) +
  labs(title="PTSD module red") + theme_classic()
df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$skyblue)
ggplot(df, aes(x=PrimaryDx, y=Eigengene)) +
  geom_boxplot() +
  # scale_fill_brewer(palette="Blues") +
  # scale_fill_manual(values=c("grey60","grey80")) +
  geom_jitter(aes(fill=PrimaryDx), color="black", size=4, shape=21, position=position_jitter(0.1)) +
  labs(title="PTSD module skyblue") + theme_classic()
# df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$royalblue)
# ggplot(df, aes(x=PrimaryDx, y=Eigengene, fill=PrimaryDx)) +
#   geom_boxplot() +
#   # scale_fill_brewer(palette="Blues") +
#   scale_fill_manual(values=c("royalblue","lightblue3")) +
#   geom_jitter(shape=16, position=position_jitter(0.1)) +
#   labs(title="PTSD module royalblue") + theme_classic()
# df <- data.frame(PrimaryDx=datMeta$PrimaryDx,Eigengene=MEs$tan)
# ggplot(df, aes(x=PrimaryDx, y=Eigengene, fill=PrimaryDx)) +
#   geom_boxplot() +
#   # scale_fill_brewer(palette="Blues") +
#   scale_fill_manual(values=c("tan","tan2")) +
#   geom_jitter(shape=16, position=position_jitter(0.1)) +
#   labs(title="PTSD module tan") + theme_classic()
dev.off()

#### protein module-CSEA barplots ####
rm(list=ls()); dev.off()
library(reshape2)
# pdf('result_0128/csea_dl_barplot.pdf')
pdf('manuscript/figures_v5/4c_csea.pdf', width=8,height=4)
load('result_0703/WGCNA_CM_dl.RData')
to_plot1 <- to_plot[apply(to_plot,1,min)<.05,]
df <- melt(to_plot)
df$Enrich <- -log10(df$value)
df <- df[df$value<.05,]
names(df)[1] <- "Module"
df$Var2 <- gsub("Oligodendrocyte", "Oligo", df$Var2)
ggplot(df, aes(x=Module, y=Enrich, fill=Module)) +
  geom_bar(color="grey50", stat="identity", position=position_dodge()) +
  geom_hline(yintercept=log10(20), linetype="dashed", color = "darkorange3") +
  xlab('') + ylab('-log10(P value)') + labs(title="MDD dlPFC") +
  scale_fill_manual(values=df$Module %>% sort %>% as.character) + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position="none") +
  facet_grid(Var2~.)
rm(list=ls())
load('result_0703/WGCNA_CP_dl.RData')
to_plot1 <- to_plot[apply(to_plot,1,min)<.05,]
df <- melt(to_plot)
df$Enrich <- -log10(df$value)
df <- df[df$value<.05,]
names(df)[1] <- "Module"
df$Var2 <- gsub("Oligodendrocyte", "Oligo", df$Var2)
ggplot(df, aes(x=Module, y=Enrich, fill=Module)) +
  geom_bar(color="grey20", stat="identity", position=position_dodge()) +
  geom_hline(yintercept=log10(20), linetype="dashed", color = "darkorange3") +
  xlab('') + ylab('-log10(P value)') + labs(title="PTSD dlPFC") +
  scale_fill_manual(values=df$Module %>% sort %>% as.character %>% unique) + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position="none") +
  facet_grid(Var2~., scales="free")
dev.off()


#### miRNA-target expression, PTSD ####
## pairs: mir589/CACNA2D1, NEGR1, OPCML, CNTN1
## mir6786/ATP6V0A1, LY6H
## mir4745/NEGR1, SH3GLB2, CNTFR, CD59
library(reshape2)
rm(list=ls())
load('DIA_new.RData')
target.mir589 <- c("CACNA2D1","NEGR1",'OPCML',"CNTN1")
target.mir6786 <- c("ATP6V0A1","LY6H")
target.mir4745 <- c("NEGR1","SH3GLB2","CNTFR","CD59")
targets <- list(MIR589=target.mir589, MIR6786=target.mir6786, MIR4745=target.mir4745)
for (i in 1:2){
  target <- targets[[i]]
  sel <- which(datMeta$Region=="dlPFC")
  df <- datExpr[match(target, datProbes$GN),sel]; datMeta1 <- datMeta[sel,]
  rownames(df) <- target
  df <- df - rowMeans(df)
  df.expr <- melt(df)
  names(df.expr) <- c("Protein","Sample","Expr")
  df.expr$Dx <- datMeta1$PrimaryDx[match(df.expr$Sample, datMeta1$Sample)]
  df.expr$miRNA <- names(targets)[i]
  if (i==1) df.comb <- df.expr else df.comb <- rbind(df.comb, df.expr) ## following
  }
### plot by miRNA
# pdf('manuscript/figures_v5/4f_targetExpr.pdf', width=7, height=3)
# ggplot(subset(df.comb, Dx!="MDD"), aes(Dx, Expr, fill=Protein)) +
#   geom_boxplot() + xlab(NULL) +
#   facet_grid(~miRNA) +
#   theme_bw()
# dev.off()

### plot genewise
pdf('manuscript/figures_v5/4f_targetExpr2.pdf', width=9, height=5)
# ggplot(subset(df.comb, Dx!="MDD"), aes(Dx, Expr, fill=Protein)) +
#   geom_boxplot() + xlab(NULL) +
#   geom_jitter(shape=16, position=position_jitter(0.1)) +
#   # facet_grid(~miRNA) +
#   facet_wrap(~miRNA+Protein) +
#   theme_bw() + theme(legend.position = "none")
ggplot(subset(df.comb, Dx!="MDD"), aes(Dx, Expr)) +
  geom_boxplot() + xlab(NULL) + ylab('Abundance') +
  geom_jitter(aes(fill=Dx), color="black", shape=21, size=4, position=position_jitter(0.1)) +
  # geom_dotplot(aes(x=Dx,y=Expr), binaxis="y", stackdir="center")+
  # geom_line(aes(x=Dx, group=Sample), alpha=0.5) +
  # facet_grid(~miRNA) +
  facet_wrap(~miRNA+Protein) +
  theme_bw()
dev.off()

#### miRNA Venn diagram ####
library(ggVennDiagram)
rm(list=ls())
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
sig01$Region <- "BA9"; sig01$Dx <- "MDD"
sig02$Region <- "BA9"; sig02$Dx <- "PTSD"
sig01.ba9 <- sig01; sig02.ba9 <- sig02
sig.ba9 <- merge(sig01, sig02[c(2,5:7)], by="Geneid")
# deg <- sig.ba9
load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
sig01$Region <- "BA25"; sig01$Dx <- "MDD"
sig02$Region <- "BA25"; sig02$Dx <- "PTSD"
sig01.ba25 <- sig01; sig02.ba25 <- sig02
sig.ba25 <- merge(sig01, sig02[c(2,5:7)], by="Geneid")
deg <- merge(sig.ba9, sig.ba25, by="Geneid")
names(deg) <- names(deg) %>% gsub('x.x','MDD.dlPFC',.) %>% gsub('y.x','PTSD.dlPFC',.) %>% gsub('x.y','MDD.sgPFC',.) %>% gsub('y.y','PTSD.sgPFC',.)
# x <- list(`PTSD dlPFC`=deg$Geneid[deg$pvalue.PTSD.dlPFC<.05] %>% as.character,
#           `PTSD sgPFC`=deg$Geneid[deg$pvalue.PTSD.sgPFC<.05] %>% as.character,
#           `MDD dlPFC`=deg$Geneid[deg$pvalue.MDD.dlPFC<.05] %>% as.character,
#           `MDD sgPFC`=deg$Geneid[deg$pvalue.MDD.sgPFC<.05] %>% as.character)
x <- list(`PTSD dlPFC`=deg$Geneid[deg$pvalue.PTSD.dlPFC<.05] %>% as.character,
          `PTSD sgPFC`=deg$Geneid[deg$pvalue.PTSD.sgPFC<.05] %>% as.character)
ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey80", high = "coral3")
pdf('manuscript/figures_v5/1_venn.pdf', width=6, height=3)
ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  scale_fill_gradient(low = "grey80", high = "coral3")
dev.off()

### style as in NN21
c3 <- deg[c("pvalue.PTSD.dlPFC","pvalue.PTSD.sgPFC")]
names(c3) <- c("PTSD dlPFC","PTSD sgPFC")
c3 <- (c3<.05)
library(limma)
a <- vennCounts(c3)
pdf('manuscript/figures_v5/1b_venn.pdf')
vennDiagram(a, circle.col = 2:3)
dev.off()

#### protein Venn diagram CMP ####
library(ggVennDiagram)
rm(list=ls())
dep.dl <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
dep.sg <- read.csv('result_0518/DE_MDD_PTSD_sg_incA_new1.csv')
dep.dl.pm <- read.csv('result_0518/DE_PvM_dl_incA_new1.csv')
dep.sg.pm <- read.csv('result_0518/DE_PvM_sg_incA_new1.csv')
dep <- merge(dep.dl, dep.dl.pm, by="EN")
x <- list(`PTSD vs CON`=dep$EN[dep$PTSD.P.Value <.05],
          `MDD vs CON`=dep$EN[dep$MDD.P.Value <.05],
          `PTSD vs MDD`=dep$EN[dep$P.Value<.05])
pdf('manuscript/supp/venn_cmp.pdf', width=12, height=9)
ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey90", high = "coral2") +
  labs(title='Venn diagram dlPFC')
dep <- merge(dep.sg, dep.sg.pm, by="EN")
x <- list(`PTSD vs CON`=dep$EN[dep$PTSD.P.Value <.05],
          `MDD vs CON`=dep$EN[dep$MDD.P.Value <.05],
          `PTSD vs MDD`=dep$EN[dep$P.Value<.05])
ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey90", high = "coral2") +
  labs(title='Venn diagram sgPFC')
dev.off()


#### miRNA Venn diagram CMP ####
library(ggVennDiagram)
rm(list=ls())
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
sig012 <- merge(sig01, sig02[c(2,5:7)], by="Geneid")
deg <- merge(sig012, sig12[c(2,5:7)], by="Geneid")
names(deg) <- names(deg) %>% gsub('.x','.MDD',.) %>% gsub('.y','.PTSD',.)
names(deg)[13:15] <- paste0(names(deg)[13:15], ".MvP")
deg$Geneid <- as.character(deg$Geneid)
x <- list(`PTSD vs CON`=deg$Geneid[deg$pvalue.PTSD<.05],
          `MDD vs CON`=deg$Geneid[deg$pvalue.MDD<.05],
          `PTSD vs MDD`=deg$Geneid[deg$pvalue.MvP<.05])
gg.dlpfc <- ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey90", high = "coral2") +
  labs(title='Venn diagram dlPFC')
load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
sig012 <- merge(sig01, sig02[c(2,5:7)], by="Geneid")
deg <- merge(sig012, sig12[c(2,5:7)], by="Geneid")
names(deg) <- names(deg) %>% gsub('.x','.MDD',.) %>% gsub('.y','.PTSD',.)
names(deg)[13:15] <- paste0(names(deg)[13:15], ".MvP")
deg$Geneid <- as.character(deg$Geneid)
x <- list(`PTSD vs CON`=deg$Geneid[deg$pvalue.PTSD<.05],
          `MDD vs CON`=deg$Geneid[deg$pvalue.MDD<.05],
          `PTSD vs MDD`=deg$Geneid[deg$pvalue.MvP<.05])
gg.sgpfc <- ggVennDiagram(x, show_intersect = F, label_alpha = 0) +
  # scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
  scale_fill_gradient(low = "grey90", high = "coral2") +
  labs(title='Venn diagram sgPFC')
pdf('manuscript/figures_v5/5c_venn_cmp_miRNA.pdf', width=8, height=6)
plot(gg.dlpfc)
plot(gg.sgpfc)
dev.off()
#### protein pathway bubble plot ####
rm(list=ls())
# mdd_dl <- xlsx::read.xlsx('result_0128/ipa/ipa_mdd_f2_canonical.xls', sheetIndex=1)
# names(mdd_dl) <- mdd_dl[1,] %>% as.matrix %>% as.character; mdd_dl <- mdd_dl[-1,]; mdd_dl$Dx <- "MDD"; mdd_dl$Region <- "dlPFC"
# ptsd_dl <- xlsx::read.xlsx('result_0128/ipa/ipa_ptsd_f2_canonical.xls', sheetIndex=1)
# names(ptsd_dl) <- ptsd_dl[1,] %>% as.matrix %>% as.character; ptsd_dl <- ptsd_dl[-1,]; ptsd_dl$Dx <- "PTSD"; ptsd_dl$Region <- "dlPFC"
# mdd_sg <- xlsx::read.xlsx('result_0128/ipa/ipa_can_mdd_sgpfc.xls', sheetIndex=1)
# names(mdd_sg) <- mdd_sg[1,] %>% as.matrix %>% as.character; mdd_sg <- mdd_sg[-1,]; mdd_sg$Dx <- "MDD"; mdd_sg$Region <- "sgPFC"
# ptsd_sg <- xlsx::read.xlsx('result_0128/ipa/ipa_can_ptsd_sgpfc.xls', sheetIndex=1)
# names(ptsd_sg) <- ptsd_sg[1,] %>% as.matrix %>% as.character; ptsd_sg <- ptsd_sg[-1,]; ptsd_sg$Dx <- "PTSD"; ptsd_sg$Region <- "sgPFC"
mdd_dl <- read.csv('result_210128/ipa/ipa_can_mdd_dl_new1.csv', skip=1)
mdd_dl$Dx <- "MDD"; mdd_dl$Region <- "dlPFC"
ptsd_dl <- read.csv('result_210128/ipa/ipa_can_ptsd_dl_new1.csv', skip=1)
ptsd_dl$Dx <- "PTSD"; ptsd_dl$Region <- "dlPFC"
mdd_sg <- read.csv('result_210128/ipa/ipa_can_mdd_sg_new1.csv', skip=1)
mdd_sg$Dx <- "MDD"; mdd_sg$Region <- "sgPFC"
ptsd_sg <- read.csv('result_210128/ipa/ipa_can_ptsd_sg_new1.csv', skip=1)
ptsd_sg$Dx <- "PTSD"; ptsd_sg$Region <- "sgPFC"

### bubble plots
df <- rbind(mdd_dl, ptsd_dl, mdd_sg, ptsd_sg)
df$`mlog(p-value)` <- df$`X.log.p.value.` %>% as.character %>% as.numeric
df$`z-score` <- df$`z.score` %>% as.character %>% as.numeric
df$Dx <- factor(df$Dx, levels=c('PTSD','MDD'))
df$Group <- paste0(df$Dx, ' ', df$Region)
shared.path <- names(table(df$`Ingenuity.Canonical.Pathways`)[table(df$`Ingenuity.Canonical.Pathways`)==4])
df <- df[df$Ingenuity.Canonical.Pathways %in% shared.path,] 
df <- df[!is.na(df$`z-score`) & df$`z-score`!=0 & df$`mlog(p-value)`>1,]
df0 <- df

pdf('manuscript/figures/5c_pathways_down.pdf', height=3, width=8.5)
## UP regulated
sel.path <- df0$Ingenuity.Canonical.Pathways[df0$Dx=="PTSD" & df0$Region=="dlPFC"][order(df0$`z-score`[df0$Dx=="PTSD"  & df0$Region=="dlPFC"], 
                                                                   decreasing = T)[1:10]]
df <- df0[df0$Ingenuity.Canonical.Pathways %in% sel.path, ]
df$Ingenuity.Canonical.Pathways <- factor(df$Ingenuity.Canonical.Pathways,
                                            ## up-regulated top 10
                                            levels=df$Ingenuity.Canonical.Pathways[df$Dx=="PTSD"&df$Region=="dlPFC"][order(df$`mlog(p-value)`[df$Dx=="PTSD"&df$Region=="dlPFC"])] %>%
                                            unique)
ggplot(df, aes(x=Dx, y=Ingenuity.Canonical.Pathways)) +
  geom_point(aes(size=`mlog(p-value)`, color=`z-score`)) +
  scale_color_gradient2(low="blue",high="red") +
  xlab('Group') + ylab('') + labs(title="IPA canonical pathway enrichment") +
  facet_grid(.~Region) + theme_bw()
dev.off()

pdf('manuscript/figures/5c_pathways_up.pdf', height=3, width=8.5)
## DOWN regulated
sel.path <- df0$Ingenuity.Canonical.Pathways[df0$Dx=="PTSD"  & df0$Region=="dlPFC"][order(df0$`z-score`[df0$Dx=="PTSD"  & df0$Region=="dlPFC"])[1:10]]
df <- df0[df0$Ingenuity.Canonical.Pathways %in% sel.path, ]
df$Ingenuity.Canonical.Pathways <- factor(df$Ingenuity.Canonical.Pathways,
                                            ## down-regulated top 10
                                            levels=df$Ingenuity.Canonical.Pathways[df$Dx=="PTSD" &df$Region=="dlPFC"][order(df$`mlog(p-value)`[df$Dx=="PTSD" &df$Region=="dlPFC"])] %>% unique)
ggplot(df, aes(x=Dx, y=Ingenuity.Canonical.Pathways)) +
  geom_point(aes(size=`mlog(p-value)`, color=`z-score`)) +
  scale_color_gradient2(low="blue",high="red") +
  xlab('Group') + ylab('') + labs(title="IPA canonical pathway enrichment") +
  facet_grid(.~Region) + theme_bw()
dev.off()

#### protein PTSD module enrichment of other traits ####
#### protein vs miRNA: SLC32A1 and miR589 ####
rm(list=ls())
getFKPM <- function(dat, log2=TRUE, minute=1e-3, trim=5, suffix=TRUE, lengths=NULL){
  sel <- (trim+1):dim(dat)[2]
  seldat <- dat[,sel]
  sums <- colSums(seldat)/1000000
  seldat <- t(t(seldat)/sums)
  len <- lengths/1000
  seldat <- seldat/len
  if(log2)  seldat <- log(seldat+minute, base=2)
  if(suffix) seldat <- cbind(dat[,1:trim], seldat)
  return(seldat)    
}
load('rawdata/dat_upmc_mirna_all1.5.RData') ## dlPFC & sgPFC
## normalization on selected miRNAs (option 1)
sel.sam <- datMeta$region=="dlPFC" ## CHECK POINT 1
sel.gene <- datProbes$TranscriptType=="miRNA"
datExpr <- datExpr[sel.gene,sel.sam]; datProbes <- datProbes[sel.gene,]; datMeta <- datMeta[sel.sam,]
datExpr.t <- datExpr; datProbes.t <- datProbes; datMeta.t <- datMeta
sel.gene2 <- rowSums(datExpr.t) > ncol(datExpr.t) * 0.5
datExpr.t <- datExpr.t[sel.gene2,]; datProbes.t <- datProbes.t[sel.gene2,]
datFPKM.t <- getFKPM(datExpr.t, trim=0, suffix=F, lengths=datProbes.t$Length)
rownames(datFPKM.t) <- datProbes.t$Geneid

## Dx specificity
sel <- datMeta.t$PrimaryDx != "MDD" ## CHECK POINT 2
datFPKM.t <- datFPKM.t[,sel]; datMeta.t <- datMeta.t[sel,]

### proteins
load('DIA_new.RData')
datExpr.p <- datExpr; datMeta.p <- datMeta; datProbes.p <- datProbes
sel.sam <- datMeta.p$Region=='dlPFC' ## CHECK POINT 3
# sel.sam <- datMeta.p$Region=='sgPFC' ## CHECK POINT 3
datExpr.p <- datExpr.p[,sel.sam]; datMeta.p <- datMeta.p[sel.sam,]
datFPKM.t1 <- datFPKM.t[,match(datMeta.p$BrNum,datMeta.t$BrNum)]
rownames(datExpr.p) <- datProbes.p$EN_short
cors <- cor(t(datFPKM.t1),t(datExpr.p), use = "complete.obs")
# cors %>% abs %>% apply(1,max) %>% hist(main="Max correlation for each miRNA", breaks=50)
# cors %>% abs %>% apply(1,median) %>% hist(main="Median correlation for each miRNA", breaks=50)
df <- melt(cors)
df$pvalue <- sapply(colnames(cors), function(g2){ 
  sapply(1:nrow(cors), function(g1){cor.test(datFPKM.t1[g1,],datExpr.p[g2,])$p.value})}) %>% melt %>% .$value

## combine with protein results
df$miRNA <- datProbes.t$Genename[match(df$Var1,datProbes.t$Geneid)]
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
# load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
df$PTSD.padj <- sig02$padj[match(df$Var1, sig02$Geneid)]
df$PTSD.pvalue <- sig02$pvalue[match(df$Var1, sig02$Geneid)]
df$MDD.padj <- sig01$padj[match(df$Var1, sig01$Geneid)]
df$MDD.pvalue <- sig01$pvalue[match(df$Var1, sig01$Geneid)]
# save(df, file="result_mirna/pair_sgpfc_CP_norm.RData")

# load('result_mirna/pair_dlpfc_CP.RData')
## run the section above
prot = "VIAAT"; mir = "MIR589"
# plot(datFPKM.t1[which(datProbes.t$Genename==mir),se[] l.sam],
#      datExpr.p[datProbes.p$EN_short==prot,sel.sam],
#      xlab="log2FPKM", ylab="log10LFQ", pch=16, col=datMeta.p$PrimaryDx[sel.sam])
df.plot <- data.frame(miRNA=datFPKM.t1[which(datProbes.t$Genename==mir),],
                      protein=datExpr.p[datProbes.p$EN_short==prot,],
                      Dx = datMeta.p$PrimaryDx)
cor(df.plot$miRNA,df.plot$protein, use = "complete.obs")
lm(miRNA~protein, data=df.plot) %>% summary
pdf('manuscript/supplementary/suppfig_10_example2.pdf', width=5, height=5)
ggplot(df.plot, aes(miRNA, protein, col="black")) +
  geom_point() + xlab('miRNA (log2FPKM)') + ylab('protein (log10LFQ)') + 
  labs(title=paste0("miRNA ", mir, " vs protein ", prot, ", cor=-0.27, R2=0.037")) +
  geom_smooth(method=lm, se=TRUE, fullrange=T, level=0.95) +
  theme_bw() + theme(legend.position="none")
dev.off()
#### miRNA-protein pairs ####
rm(list=ls())
##MDD
pairs <- read.csv('manuscript/results/pair_dlpfc_mdd_share_511_enrichment_new1.csv')
mirna_sig <- pairs$miRNA[pairs$p.chisq<.05]
modules <- c("darkred","grey60")
load('result_mirna/pair_dlpfc_CM.RData')
df.mod <- subset(df, PTSD.mod %in% modules & miRNA %in% mirna_sig & value < 0 & pvalue < 0.05 & MDD.pvalue < 0.05)
dim(df.mod)
load('DIA_new.RData')
df.mod$Genename <- datProbes$GN[match(df.mod$Var2,datProbes$EN_short)]
names(df.mod) <- c("Geneid","Protein","Cor","Pvalue.Cor","miRNA","miRNA.Pvalue.MDD","miRNA.Pvalue.PTSD","Module","Genename")
df.mod$PrimaryDx <- "MDD"
df.mod.mdd <- df.mod
rm(datExpr,datMeta,datProbes,df,df.mod,pairs,mirna_sig,modules)

##PTSD
pairs <- read.csv('manuscript/results/pair_dlpfc_ptsd_share_511_enrichment_new1.csv')
mirna_sig <- pairs$miRNA[pairs$p.chisq<.05]
modules <- c("skyblue","red")
load('result_mirna/pair_dlpfc_CP.RData')
df.mod <- subset(df, PTSD.mod %in% modules & miRNA %in% mirna_sig & value < 0 & pvalue < 0.05 & PTSD.pvalue < 0.05)
dim(df.mod)
load('DIA_new.RData')
df.mod$Genename <- datProbes$GN[match(df.mod$Var2,datProbes$EN_short)]
names(df.mod) <- c("Geneid","Protein","Cor","Pvalue.Cor","miRNA","miRNA.Pvalue.MDD","miRNA.Pvalue.PTSD","Module","Genename")
df.mod$PrimaryDx <- "PTSD"
df.mod.ptsd <- df.mod
rm(datExpr,datMeta,datProbes,df,df.mod,pairs,mirna_sig,modules)
df.mod <- rbind(df.mod.ptsd, df.mod.mdd)
df.mod <- df.mod[c(10,8,1,9,2,5,3,4,6,7)]
write.csv(df.mod, row.names=F, file="manuscript/results/pairs_sig_mirna_protein_2dx2mod.csv")

##for mRNA-protein divergent genes
rm(list=ls())
pairs <- read.csv('manuscript/results/pair_dlpfc_ptsd_share_511_enrichment_new1.csv')
mirna_sig <- pairs$miRNA[pairs$p.chisq<.05]
genes <- c("RHOB","PTPRE","DPP6","CD59")
load('result_mirna/pair_dlpfc_CP.RData')
df.mod <- subset(df, Var2 %in% genes & value < 0 & pvalue < 0.05)
dim(df.mod)
load('DIA_new.RData')
df.mod$Genename <- datProbes$GN[match(df.mod$Var2,datProbes$EN_short)]
names(df.mod) <- c("Geneid","Protein","Cor","Pvalue.Cor","miRNA","miRNA.Pvalue.MDD","miRNA.Pvalue.PTSD","Module","Genename")
df.mod$PrimaryDx <- "PTSD"
df.mod <- df.mod[c(10,8,1,9,2,5,3,4,6,7)]
write.csv(df.mod, file="manuscript/results/pairs_divergent-proteins.csv", row.names=F)





#### histogram: RNA biotype fractions ####
rm(list=ls())
getFKPM <- function(dat, log2=TRUE, minute=1e-3, trim=5, suffix=TRUE){
  sel <- (trim+1):dim(dat)[2]
  seldat <- dat[,sel]
  sums <- colSums(seldat)/1000000
  seldat <- t(t(seldat)/sums)
  len <- dat$Length/1000
  seldat <- seldat/len
  if(log2)  seldat <- log(seldat+minute, base=2)
  if(suffix) seldat <- cbind(dat[,1:trim], seldat)
  return(seldat)    
}
# load('../data/dat_upmc_r_keep.RData')
load('../transcriptomics_ctx/data/dat_upmc_r.RData')
fpkm_new <- getFKPM(dat13r, trim=6, suffix=F)
fpkm_rm <- rowMeans(fpkm_new)
hist(fpkm_rm[fpkm_rm>log2(1e-3)], breaks=40)
hist(fpkm_rm[fpkm_rm>log2(1e-3) & suff$Biotype=="protein_coding"], breaks=40)
load('DIA_new.RData')
gene_prot <- datProbes$GN[datProbes$DB=="sp"]
hist(fpkm_rm[fpkm_rm>log2(1e-3) & suff$Genename %in% gene_prot], breaks=40)
df <- data.frame(Genename=suff$Genename, Type="All transcript", FPKM=fpkm_rm)
df$Type <- as.character(df$Type)
df$Type[which(suff$Biotype=="protein_coding")] <- "Coding transcript"
df$Type[df$Genename %in% gene_prot] <- "Protein"
df$Color <- df$Type %>% gsub("All","grey",.) %>% gsub("Coding","black",.) %>% gsub("Protein","darkred",.) %>% gsub(" transcript","",.)
gene_exp <- rowSums(dat13r[-(1:6)]>0)>ncol(dat13r[-(1:6)])/100 ##5% presence
fill_color <- df$Color %>% as.character
# ggplot(df[gene_exp,], aes(x=FPKM, fill=Type)) +
#   geom_histogram(color="white", alpha=0.9, position="stack", bins=80) +
#   theme_classic()
# dev.off()
library(ggplot2)
pdf('manuscript/figures/1_RNAtype.pdf', height=5, width=8)
ggplot(df[gene_exp,], aes(x=FPKM, fill=Type)) +
  geom_histogram(color="white",alpha=0.9, position="stack", binwidth = 0.3) +
  scale_fill_discrete(type = c('grey80','grey20','darkred')) + xlab('log2(FPKM)') + ylab('Counts') +
  theme(text = element_text(size=10), legend.position = "right") + xlim(c(-10,10)) +
  theme_classic()
dev.off()
## Tech adjustments>>>
## capped at 10
## using UPMC transcriptomic reference
## FPKM grounded at 1e-3
## sp proteins only


#### RRHO of Dx x region ####
library(RRHO) ## sample code below
# list.length <- 100
# list.names <- paste('Gene',1:list.length, sep='')
# gene.list1<- data.frame(list.names, sample(-100:100, 100))
# gene.list2<- data.frame(list.names, sample(-100:100, 100))
# # Enrichment alternative
# RRHO.example <-  RRHO(gene.list1, gene.list2, alternative='enrichment')
# image(RRHO.example$hypermat)
# 
# # Two sided alternative
# RRHO.example <-  RRHO(gene.list1, gene.list2, alternative='two.sided')
# image(RRHO.example$hypermat)

### 
# dep.dl <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
# dep.sg <- read.csv('result_0518/DE_MDD_PTSD_sg.csv')
dep.dl <- read.csv('result_200518/DE_MDD_PTSD_dl_incA_new1.csv')
dep.sg <- read.csv('result_200518/DE_MDD_PTSD_sg_incA_new1.csv')
dep <- merge(dep.dl, dep.sg, by="EN")
names(dep) <- names(dep) %>% gsub('x','dl',.) %>% gsub('y','sg',.)
write.csv(dep[-c(2,21)], file="result_200518/DE_MDD_PTSD_dl_sg.csv", row.names=F)
dep$rank.mdd.dl <- rank(dep$MDD.logFC.dl)
dep$rank.ptsd.dl <- rank(dep$PTSD.logFC.dl)
dep$rank.mdd.sg <- rank(dep$MDD.logFC.sg)
dep$rank.ptsd.sg <- rank(dep$PTSD.logFC.sg)
# dep.sel <- dep[c(1,20,40,41,3,11)]
# write.table(dep.sel, file="result_0518/DE_rrho_dl.txt", row.names = F, quote=F)
pdf('manuscript/figures/1_rrho.pdf')
RRHO.example <-  RRHO(dep[c(1,40)], dep[c(1,41)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="MDD dlPFC", ylab="PTSD dlPFC")
RRHO.example <-  RRHO(dep[c(1,42)], dep[c(1,41)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="MDD sgPFC", ylab="PTSD dlPFC")
RRHO.example <-  RRHO(dep[c(1,43)], dep[c(1,41)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="PTSD sgPFC", ylab="PTSD dlPFC")
RRHO.example <-  RRHO(dep[c(1,42)], dep[c(1,40)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="MDD sgPFC", ylab="MDD dlPFC")
RRHO.example <-  RRHO(dep[c(1,43)], dep[c(1,40)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="PTSD sgPFC", ylab="MDD dlPFC")
RRHO.example <-  RRHO(dep[c(1,42)], dep[c(1,43)], alternative='two.sided')
image(RRHO.example$hypermat, col=hcl.colors(12, "RdYlBu", rev = T), 
      xlab="MDD sgPFC", ylab="PTSD sgPFC")
dev.off()


#### MDD/PTSD DEP CSEA in barplots [21/9/14] ####
rm(list=ls())
# BiocManager::install('gdata')
# install.packages(type="source", repos=NULL, "~/Downloads/pSI_1.1.tar.gz")
library(pSI)
library(gplots)
library(ggplot2)
library(WGCNA)
# sig <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
# sig <- read.csv('result_0518/DE_MDD_PTSD_sg.csv')
sig <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
# sig <- read.csv('result_0518/DE_MDD_PTSD_sg_incA_new1.csv')
load('../data/datProbes.RData')
sig$Geneid <- datProbes$ensembl_gene_id[match(sig$GN, datProbes$external_gene_id)]
load('../data/zhang.pSIout.RData')
gene.uniq <- intersect(sig$Geneid, rownames(pSI.output)) %>% unique
pSI.output <- pSI.output[rownames(pSI.output) %in% gene.uniq,]
# sig <- sig[sig$Geneid %in% gene.uniq & !is.na(sig$PTSD.P.Value),] ## PTSD
sig <- sig[sig$Geneid %in% gene.uniq & !is.na(sig$MDD.P.Value),] ## MDD
# cats <- paste0('Top',1:10,'00')
cats <- paste0('Pvalue<',0.05*1:10)
ncat <- length(cats)
cell.p.zhang = matrix(NA, ncat, 5);  rownames(cell.p.zhang) = cats
colnames(cell.p.zhang) = colnames(pSI.output)
for(j in 1:2) {
  ## ordered by top ith 100 genes
  # f = fisher.iteration(pSI.output, sig$Geneid[order(sig$PTSD.P.Value, decreasing = F)][(100*(j-1)+1):(100*j)], p.adjust = F)
  # f = fisher.iteration(pSI.output, sig$Geneid[order(sig$PTSD.P.Value, decreasing = F)][1:(100*j)], p.adjust = F)
  ## ordered by P values
  ## PTSD
  # f = fisher.iteration(pSI.output, sig$Geneid[sig$PTSD.P.Value < 0.05 * j & sig$PTSD.P.Value >= 0.05 * (j-1)], p.adjust = F)
  ## MDD
  f = fisher.iteration(pSI.output, sig$Geneid[sig$MDD.P.Value < 0.05 * j & sig$MDD.P.Value >= 0.05 * (j-1)], p.adjust = F)
  cell.p.zhang[j,] = f$`0.05 - nominal`
}
cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
to_plot = cell.p.zhang.fdr %>% as.data.frame
dendro.col = dendro.col.zhang
to_plot$pvalue <- rownames(to_plot)
df <- reshape2::melt(to_plot)
df$Dx <- "MDD"; df$BA <- "dlPFC"
write.csv(df, row.names=F, file="result_0128/csea_cm_dl_p05_incA.csv")

## plot
rm(list=ls())
# df.cm.sg <- read.csv('result_0128/csea_cm_sg_p05.csv')
# df.cm.dl <- read.csv('result_0128/csea_cm_dl_p05.csv')
# df.cp.sg <- read.csv('result_0128/csea_cp_sg_p05.csv')
# df.cp.dl <- read.csv('result_0128/csea_cp_dl_p05.csv')
df.cm.sg <- read.csv('result_0128/csea_cm_sg_p05_incA.csv')
df.cm.dl <- read.csv('result_0128/csea_cm_dl_p05_incA.csv')
df.cp.sg <- read.csv('result_0128/csea_cp_sg_p05_incA.csv')
df.cp.dl <- read.csv('result_0128/csea_cp_dl_p05_incA.csv')
df <- rbind(df.cm.sg, df.cm.dl, df.cp.sg, df.cp.dl)
df$pvalue <- df$pvalue %>% gsub("value","",.) %>% gsub("<"," < ",.)
names(df)[1] <- "P-value"
pdf('manuscript/figures/1_depCSEA.pdf', width=9)
ggplot(df[df$`P-value` %in% c("P < 0.05","P < 0.1"),], 
       aes(variable, -log10(value), fill=`P-value`)) + 
  geom_bar(stat="identity", position="dodge", color="black") +
  geom_hline(yintercept=-log10(0.5), col="grey40", lty=2) +
  facet_grid(rows = vars(BA), cols = vars(Dx), switch = "y") +
  xlab(NULL) + ylab('-log10(P-value)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8)) 
dev.off()




#### how smRNA differs from bulk RNAseq ####
library(reshape2)
rm(list=ls())
getFKPM <- function(dat, log2=TRUE, minute=1e-3, trim=5, suffix=TRUE, lengths=NULL){
  sel <- (trim+1):dim(dat)[2]
  seldat <- dat[,sel]
  sums <- colSums(seldat)/1000000
  seldat <- t(t(seldat)/sums)
  len <- lengths/1000
  seldat <- seldat/len
  if(log2)  seldat <- log(seldat+minute, base=2)
  if(suffix) seldat <- cbind(dat[,1:trim], seldat)
  return(seldat)    
}
load('rawdata/dat_upmc_mirna_all1.5.RData') ## dlPFC & sgPFC
## normalization on selected miRNAs (option 1)
# sel.sam <- datMeta$region=="dlPFC" &## CHECK POINT 1
sel.sam <- 1:nrow(datMeta)
sel.gene <- datProbes$TranscriptType=="miRNA"
datExpr <- datExpr[sel.gene,sel.sam]; datProbes <- datProbes[sel.gene,]; datMeta <- datMeta[sel.sam,]
datExpr.t <- datExpr; datProbes.t <- datProbes; datMeta.t <- datMeta
sel.gene2 <- rowSums(datExpr.t) > ncol(datExpr.t) * 0.5
datExpr.t <- datExpr.t[sel.gene2,]; datProbes.t <- datProbes.t[sel.gene2,]
datFPKM.t <- getFKPM(datExpr.t, trim=0, suffix=F, lengths=datProbes.t$Length)
rownames(datFPKM.t) <- datProbes.t$Geneid
load('../data/dat_upmc_r.RData')
shared <- intersect(datProbes.t$Geneid, suff$Geneid)
id1 <- match(shared, datProbes.t$Geneid); datFPKM.t <- datFPKM.t[id1,]; datProbes.t <- datProbes.t[id1,]
id2 <- match(shared, suff$Geneid); fpkm.m <- fpkm[id2,[]]; suff.m <- suff[id2,]
shared.sam <- match(datMeta.t$Sample, meta13r$Sample); fpkm.m <- fpkm.m[,shared.sam]; datMeta.m <- meta13r[shared.sam,]
df.plot <- data.frame(bulk=rowMeans(fpkm.m), smRNA=rowMeans(datFPKM.t))
summary(lm(smRNA~bulk, df.plot))
pdf('manuscript/supp/bulk_vs_smRNA.pdf')
ggplot(df.plot, aes(bulk, smRNA)) +
  geom_point() + 
  xlab('mean log2FPKM (bulk RNA-seq)') + ylab('mean log2FPKM (smRNA-seq)') + 
  labs(title="Bulk vs smRNA-seq, R2=0.0063") +
  geom_smooth(method=lm, se=TRUE, fullrange=T, level=0.95) +
  theme_bw()
dev.off()






#### enrichment: DE miRNA in protein/mRNA module ####
rm(list=ls())
# load('result_mirna/pair_dlpfc_CP.RData') ## CP dlPFC
# load('result_mirna/pair_dlpfc_CM.RData') ## CM dlPFC
load('result_mirna/pair_sgpfc_CP_norm.RData') ## CP sgPFC
# load('result_mirna/pair_sgpfc_CM_norm.RData') ## CM sgPFC
## test module-wise enrichment
dim(df) ## from pair-wise correlation results
# load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData') ## BA9
load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData') ## BA25
df$MDD.pvalue <- sig01$pvalue[match(df$Var1,sig01$Geneid)]
df$PTSD.pvalue <- sig02$pvalue[match(df$Var1,sig02$Geneid)]
## swtich between PTSD and MDD
# load('result_0703/WGCNA_CP_dl.RData') ## PTSD dlPFC
# load('result_0703/WGCNA_CM_dl.RData') ## MDD dlPFC
load('result_0703/WGCNA_CP_sg.RData') ## PTSD sgPFC
# load('result_0703/WGCNA_CM_sg.RData') ## MDD sgPFC
# all <- all[all$Protein %in% datProbes.p$EN,]
all$Protein <- all$Protein %>% gsub("_HUMAN","",.)
df$PTSD.mod <- all$module[match(df$Var2,all$Protein)] ## name PTSD.mod doesn't matter
sel.m <- df$PTSD.pvalue<.01 & df$pvalue<.05 ## PTSD
# sel.m <- df$MDD.pvalue<.01 & df$pvalue<.01 ## MDD
df1 <- df[which(sel.m),]
dim(df1)

# ## calculating odds ratio
# or0 <- table(df$PTSD.mod) %>% as.vector(.)/nrow(df)
# freq <- df1 %>% count(Var1, PTSD.mod) %>% tidyr::spread(key="Var1",value="n")
# rownames(freq) <- freq$PTSD.mod
# freq <- freq[-1]; freq[is.na(freq)] <- 0
# freq <- freq/colSums(freq)
# odds <- freq/or0
# odds$module <- rownames(odds)
# df.p <- melt(odds)
# df.p$Genename <- datProbes.t$Genename[match(df.p$variable,datProbes.t$Geneid)]
# df.p$Genename[is.na(df.p$Genename)] <- df.p$variable[is.na(df.p$Genename)]
# df.p$score <- ""
# df.p$score[df.p$value>3] <- ">3"
# df.p$score[df.p$value>5] <- ">5"
# df.p$score[df.p$value>10] <- ">10"

## calculating OR chisq test
freq <- df1 %>% count(Var1, PTSD.mod) %>% tidyr::spread(key="Var1",value="n")
rownames(freq) <- freq$PTSD.mod
freq <- freq[-1]; freq[is.na(freq)] <- 0
freq$module <- rownames(freq)
df.p <- reshape2::melt(freq)
# df.p$nmod <- table(df1$PTSD.mod) %>% as.vector
# df.p$nrna <- table(factor(df1$Var1)) %>% as.vector %>% rep(each=length(table(df1$PTSD.mod)))
# df.p$exp <- df.p$nrna * (table(df1$PTSD.mod) %>% as.vector)/nrow(df1)
df.p$n22 <- sum(df.p$value) - df.p$value
df.p$n12 <- (table(factor(df1$Var1)) %>% as.vector %>% rep(each=length(table(factor(df1$PTSD.mod))))) -df.p$value
df.p$n21 <- (table(factor(df1$PTSD.mod)) %>% as.vector) - df.p$value
df.p$p.chisq <- apply(df.p[c(3,5,6,4)], 1, function(d) chisq.test(matrix(d %>% as.vector, ncol=2))$p.value)
df.p$p.sign <- sign(df.p$value*df.p$n22/(df.p$n12*df.p$n21)-1)
# df.p$Genename <- datProbes.t$Genename[match(df.p$variable,datProbes.t$Geneid)]
# df.p$Genename[is.na(df.p$Genename)] <- df.p$variable[is.na(df.p$Genename)]
df.p$Genename <- sig01$Transcriptid[match(df.p$variable,sig01$Geneid)] %>% as.character
df.p$Genename[df.p$Genename==""] <- df.p$variable[df.p$Genename==""] %>% as.character
df.p$score <- ""
df.p$score[df.p$p.chisq<.05 & df.p$p.sign > 0] <- "*"
df.p$score[df.p$p.chisq<.01 & df.p$p.sign > 0] <- "**"
df.p$score[df.p$p.chisq<.001 & df.p$p.sign > 0] <- "***"


## plot heatmap
library(reshape2)
library(ggplot2)
pdf('manuscript/supp/enrichment_mirna_module_sg_cp.pdf', width=8, height=3)
ggplot(df.p[df.p$module!="grey",], aes(Genename, module, fill=value)) + 
  geom_tile() +
  geom_text(aes(label = score), vjust = .8) + 
  scale_fill_gradient2(low="blue", mid="white", high="firebrick2") +
  xlab('miRNAs') + ylab('Protein modules') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=8))
dev.off()

## quick check
sig01[sig01$Transcriptid %in% c("MIR218-1","MIR379","MIR589","MIR181A1","MIR107"),] %>% View()
sig02[sig02$Transcriptid %in% c("MIR218-1","MIR379","MIR589","MIR199A2"),] %>% View()




#####