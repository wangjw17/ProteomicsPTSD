#### Proteomics analysis, Apr 1, 2020 

### read in data ====
dat <- xlsx::read.xlsx('Girgenti merged analysis/Samples Table of Girgenti Human Exp1 and Exp2 log10.xlsx', sheetIndex=1)
labels <- dat[2,]
dat <- dat[-(1:2),]
names(dat) <- labels %>% as.matrix %>% as.vector
View(dat)
samples <- names(dat)[-(1:8)]
samples <- gsub("_Girgenti|.mzML|_DIA","",samples)
samples <- gsub("OTF18", "sgPFC", samples)
samples <- gsub("OTF19", "dlPFC", samples)
samples <- gsub("-", "_", samples)
names(dat)[-(1:8)] <- samples
ntrim <- 8

library(stringr)
datProbes <- dat[1:ntrim]
datProbes$DB <- datProbes$`Protein Name` %>% as.character %>% strsplit(split = "\\|") %>% sapply("[",1)
datProbes$DB_Group <- 0; datProbes$DB_Group[datProbes$DB %>% grep("Group",.)] <- 1
datProbes$DB_DECOY <- 0; datProbes$DB_DECOY[datProbes$DB %>% grep("DECOY",.)] <- 1
datProbes$DB <- datProbes$DB %>% gsub("DECOY_|Group of ","",.)
datProbes$UI <- datProbes$`Protein Name` %>% as.character %>% strsplit(split = "\\|") %>% sapply("[",2)
datProbes$EN <- datProbes$`Protein Name` %>% as.character %>% strsplit(split = "\\|") %>% 
  sapply("[",3) %>% strsplit(split = " ") %>% sapply("[",1)
datProbes$EN_short <- datProbes$EN %>% gsub("_HUMAN","",.)
datProbes$PN <- datProbes$`Protein Name` %>% as.character %>% strsplit(split = "HUMAN ") %>% 
  sapply("[",2) %>% strsplit(split=" OS=") %>% sapply("[",1)
datProbes$SV = datProbes$`Protein Name` %>% str_extract(pattern = "SV=[0-9]+") %>% gsub("SV=","",.)
datProbes$PE = datProbes$`Protein Name` %>% str_extract(pattern = "PE=[0-9]+") %>% gsub("PE=","",.)
datProbes$GN = datProbes$`Protein Name` %>% str_extract(pattern = "GN=\\S+") %>% gsub("GN=","",.)
datProbes$OX = datProbes$`Protein Name` %>% str_extract(pattern = "OX=[0-9]+") %>% gsub("OX=","",.)
datProbes$OS = datProbes$`Protein Name` %>% str_extract(pattern = "OS=.+ OX") %>% gsub("OS=| OX","",.)

datMeta <- data.frame(Sample.raw=samples)
datMeta$Sample <- datMeta$Sample.raw %>% as.character %>% strsplit(split="_") %>% sapply("[",3)
datMeta$PrimaryDx <- gsub("[0-9]+","",datMeta$Sample)
datMeta$Region <- datMeta$Sample.raw %>% as.character %>% strsplit(split="_") %>% sapply("[",1)

datExpr <- dat[-(1:ntrim)]
save(datExpr, datProbes, datMeta, file="DIA_new0.RData")

## find matched information
# load("DIA_new.RData")
datExpr1 <- datExpr; 
datMeta1 <- datMeta; datProbes1 <- datProbes
load("DIA_dl_sgPFC.RData")
label <- datMeta$Sample.raw
label <- gsub("_Girgenti|_DIA|.mzML","",label)
label <- gsub("OTF18-","sgPFC_",label)
label <- gsub("OTF19-","dlPFC_",label)
datMeta$Sample.raw <- label
datMeta2 <- datMeta[match(datMeta1$Sample.raw,datMeta$Sample.raw),]
sam.id <- !is.na(datMeta2$Sample.raw)
datMeta1 <- datMeta2[sam.id,]
datExpr1 <- datExpr1[,sam.id]
datExpr <- datExpr1; datMeta <- datMeta1; datProbes <- datProbes1
save(datExpr, datProbes, datMeta, file="DIA_new.RData")


### data inspection ====
rm(list=ls())
load('DIA_new.RData')
datExpr <- datExpr %>% apply(., 2, function(x) x %>% as.character %>% as.numeric)
sum(is.na(apply(datExpr,1,sum))); sum(is.na(datExpr))
sel.rm <- is.na(apply(datExpr,1,sum))
datExpr <- datExpr[!sel.rm,]; datProbes <- datProbes[!sel.rm,]
# datMeta$Region <- datMeta$Sample.raw %>% as.character %>% strsplit(split="_") %>% sapply("[",1)
var.sel <- apply(datExpr,1,sd) %>% order(decreasing = T)
# colnames(datExpr) <- paste0(colnames(datExpr), "_", datMeta$Region)

pdf('result_0518/result_dist.pdf')
heatmap(datEpr[var.sel[1:1000],], main="Heatmap of top 1k most variated proteins")
plot(density(rowMeans(datExpr)), main="Distribution of mean expression levels")
lines(density(rowMeans(datExpr[,datMeta$Region=="dlPFC"])), col="red")
lines(density(rowMeans(datExpr[,datMeta$Region=="sgPFC"])), col="blue")
legend('topleft', legend = c("Mixed","dlPFC","sgPFC"), lty=1, col=c("black","red","blue"))
plot(density(apply(datExpr,1,var) %>% log), main="Distribution of log variance")
plot(apply(datExpr,1,mean),apply(datExpr,1,sd), main="Mean vs Sd")
dev.off()

# id <- datMeta$Region=="dlPFC"
# pca.tpm <- prcomp(t(datExpr[,id]), scale. = T, center = T)
pca.tpm <- prcomp(t(datExpr), scale. = T, center = T)
pcatpm <- pca.tpm$x
pcatsne <- tsne::tsne(pcatpm[,1:10])

pdf('result_0518/result_dimred.pdf', width=10)
par(mfrow=c(2,3))
plot(pcatpm[,1:2], main="PCA of samples", pch=16, col=factor(datMeta$Region))
legend('topleft', levels(factor(datMeta$Region)), pch=16, col=1:2, cex=.5)
plot(pcatpm[,2:3], main="PCA of samples", pch=16, col=factor(datMeta$PrimaryDx))
legend('topleft', levels(factor(datMeta$PrimaryDx)), pch=16, col=1:3, cex=.5)
plot(pcatpm[,2:3], main="PCA of samples", pch=16, col=factor(datMeta$Sex))
legend('topleft', levels(factor(datMeta$Sex)), pch=16, col=1:2, cex=.5)
plot(pcatpm[,2:3], main="PCA of samples", pch=16, col=factor(datMeta$Race))
legend('topleft', levels(factor(datMeta$Race)), pch=16, col=1:3, cex=.5)
plot(pcatpm[,2:3], main="PCA of samples", pch=16, col=factor(datMeta$Smoking))
legend('topleft', levels(factor(datMeta$Smoking)), pch=16, col=1:3, cex=.5)
plot(pcatpm[,2:3],  main="PCA of samples", pch=16, col=factor(datMeta$BrNum), )
legend('topleft', levels(factor(datMeta$BrNum)), pch=16, col=1:57, cex=.5)
par(mfrow=c(1,1))
plot(pca.tpm$sdev^2/sum(pca.tpm$sdev^2))
par(mfrow=c(2,3))
plot(pcatsne, col=factor(datMeta$Region), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$Region)), pch=16, col=1:2, cex=.5)
plot(pcatsne, col=factor(datMeta$PrimaryDx), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$PrimaryDx)), pch=16, col=1:3, cex=.5)
plot(pcatsne, col=factor(datMeta$Race), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$Race)), pch=16, col=1:3, cex=.5)
plot(pcatsne, col=factor(datMeta$Sex), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$Sex)), pch=16, col=1:2, cex=.5)
plot(pcatsne, col=factor(datMeta$Smoking), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$Smoking)), pch=16, col=1:3, cex=.5)
plot(pcatsne, col=factor(datMeta$BrNum), pch=16, main="tSNE of samples")
legend('topleft', levels(factor(datMeta$BrNum)), pch=16, col=1:57, cex=.5)
dev.off()

pdf('result_0518/result_hclust.pdf', width=12, height=8)
d <- dist(t(datExpr[var.sel,]), method = "euclidean")
hc1 <- hclust(d, method = "average")
plot(hc1, cex = 0.4, hang = -1, main="All proteins")
d <- dist(t(datExpr[var.sel[1:1000],]), method = "euclidean")
hc1 <- hclust(d, method = "average" )
plot(hc1, cex = 0.4, hang = -1, main="Top 1k proteins")
d <- dist(t(datExpr[var.sel[1:500],]), method = "euclidean")
hc1 <- hclust(d, method = "average" )
plot(hc1, cex = 0.4, hang = -1, main="Top 500 proteins")
dev.off()

library(RColorBrewer)
library(gplots)
demos <- datMeta[c('AgeDeath','PMI','RIN1')]
demos$AgeDeath <- as.numeric(demos$AgeDeath)
demos$PMI <- as.numeric(demos$PMI)
demos$RIN1 <- as.numeric(demos$RIN1)
demos$MDD <- as.numeric(datMeta$PrimaryDx=="MDD"); demos$MDD[datMeta$PrimaryDx=="PTSD"] <- NA
demos$PTSD <- as.numeric(datMeta$PrimaryDx=="PTSD"); demos$PTSD[datMeta$PrimaryDx=="MDD"] <- NA
demos$Race <- as.numeric(datMeta$Race=="W"); demos$Race[datMeta$Race=="Asian"] <- NA
demos$Sex <- as.numeric(datMeta$Sex=="F")
demos$Region <- as.numeric(datMeta$Region=="dlPFC")
pdf('result_cor_traits.pdf')
corrs <- cor(pcatpm[,1:10], demos[c('Sex','MDD','PTSD','AgeDeath','PMI','RIN1','Race','Region')], 
             use="pairwise.complete.obs")
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(abs(corrs), dendrogram = 'none', trace="none", Rowv = NULL, Colv=NULL, 
          density.info="none", col=brewer.pal(9,'YlOrRd'), cexCol = .9, cexRow = .9)
dev.off()

### differential expression analysis ====
load('DIA_new.RData')
datExpr <- datExpr %>% apply(., 2, function(x) x %>% as.character %>% as.numeric)
sum(is.na(apply(datExpr,1,sum))); sum(is.na(datExpr))
sel.rm <- is.na(apply(datExpr,1,sum))
datExpr <- datExpr[!sel.rm,]; datProbes <- datProbes[!sel.rm,]
datMeta$PrimaryDx <- datMeta$PrimaryDx %>% as.factor
datMeta$Region <- datMeta$Region %>% as.factor
datMeta$AgeDeath <- datMeta$AgeDeath %>% as.numeric
datMeta$PMI <- datMeta$PMI %>% as.numeric
save(datExpr, datMeta, datProbes, file="DIA_new.RData")

sel0 <- datMeta$Race == "Asian"
datMeta <- datMeta[!sel0,]; datExpr <- datExpr[,!sel0]
datMeta$Race <- datMeta$Race %>% factor

library(limma)
sumstats = vector(mode="list", length=2); names(sumstats)= c("MDD", "PTSD")
sam.sel <- which(datMeta$Region=="dlPFC")
mod = model.matrix(~PrimaryDx+AgeDeath+Race+Sex, data=datMeta[sam.sel,])
fit = eBayes(lmFit(datExpr[,sam.sel], mod),trend = T, robust=T)
sumstats$MDD = topTable(fit, coef=2, number = Inf, sort.by = "none",confint = T)
sumstats$PTSD = topTable(fit, coef=3, number = Inf, sort.by = "none", confint = T)
prot1 = do.call("cbind", sumstats)
prot1$PN <- datProbes$PN; prot1$GN <- datProbes$GN; prot1$EN<- datProbes$EN
View(prot1)
sum(prot1$MDD.adj.P.Val<.05); sum(prot1$PTSD.adj.P.Val<.05)
write.csv(file="result_0518/DE_MDD_PTSD_dl.csv", prot1)


## ComBat corrected
sumstats = vector(mode="list", length=2); names(sumstats)= c("MDD", "PTSD")
sam.sel <- 1:nrow(datMeta)
mod = model.matrix(~PrimaryDx+AgeDeath+Race+Sex, data=datMeta[sam.sel,])
fit = eBayes(lmFit(datExpr[,sam.sel], mod),trend = T, robust=T)
sumstats$MDD = topTable(fit, coef=2, number = Inf, sort.by = "none",confint = T)
sumstats$PTSD = topTable(fit, coef=3, number = Inf, sort.by = "none", confint = T)
prot1 = do.call("cbind", sumstats)
prot1$PN <- datProbes$PN; prot1$GN <- datProbes$GN; prot1$EN<- datProbes$EN
View(prot1)
write.csv(file="result_0518/DE_MDD_PTSD_combat.csv", prot1)


### WGCNA ====
library(magrittr)
library(stats)
library(WGCNA)
library(flashClust)
library(biomaRt)
library(corrplot)
library(ggplot2)
library(cqn)
library(sva)
library(limma)
library(statmod)
library(dplyr)
library(pSI)
library(pSI.data)
library(gplots)
library(stringr)
# options(stringsAsFactors = FALSE, digits = 3)
# theme_update(plot.title = element_text(hjust = 0.5))
rm(list=ls())
# tag = "dl" ##Approach 0: not adjusted for covariates
# tag = "dl_1" ##Approach 1: adjusted for covariates
# tag = "dl_2" ##Approach 2: not adjusted but associated calculated
# new: not adjusted for covariates

## load & prepare data
# load('DIA_dl_sgPFC.RData')
load('DIA_new.RData')

# ## select invariant proteins
# datExpr0 <- datExpr; datExpr <- datExpr[inv.sel,]
# datProbes0 <- datProbes; datProbes <- datProbes[inv.sel,]

sel0 <- !is.na(datProbes$GN)
sam0 <- datMeta$Region == "sgPFC" ## dlPFC
datExpr <- datExpr[sel0,sam0]
datMeta <- datMeta[sam0,]
datProbes <- datProbes[sel0,]
datMeta$AgeDeath <- as.numeric(datMeta$AgeDeath)
datMeta$PMI <- as.numeric(datMeta$PMI)
datMeta$Sample.disam <- paste0(datMeta$Sample,"_",datMeta$BrNum)
colnames(datExpr) <- datMeta$Sample.disam
sam1 <- datMeta$PrimaryDx!="MDD"
# sam1 <- datMeta$PrimaryDx!="PTSD"
datExpr <- datExpr[,sam1]
datMeta <- datMeta[sam1,]

##combined
# load('DIA_dl_sgPFC.RData') 
# sel0 <- !is.na(datProbes$GN)
# datExpr <- datExpr[sel0,]
# datProbes <- datProbes[sel0,]
# datMeta$Sample <- datMeta$Sample.raw %>% gsub("_","",.) %>% str_extract("CON[0-9]+|MDD[0-9]+|PTSD[0-9]+")
# names(datExpr) <- paste0(datMeta$Sample, "_",datMeta$Region)
# datMeta$Region <- factor(datMeta$Region)
# datMeta$AgeDeath <- as.numeric(datMeta$AgeDeath)
# datMeta$PMI <- as.numeric(datMeta$PMI)
# datMeta$Sample.disam <- paste0(datMeta$Sample,"_",datMeta$Region)
# colnames(datExpr) <- datMeta$Sample.disam

## remove covariates
dat <- datExpr %>% t()
datMeta0 <- datMeta[c('AgeDeath', 'RIN', 'PMI', 'Sex')]
datMeta0$Sex <- as.numeric(datMeta0$Sex)
# lm.cov <- lm(dat ~ AgeDeath + RIN1 + PMI + Sex, data=datMeta0) ##for sex specific, not including Race/Area
# dat.ori <- dat
# dat <- dat - as.matrix(datMeta0[c('AgeDeath', 'RIN1', 'PMI', 'Sex')]) %*% lm.cov$coefficients[2:5,]
##combined
# dat <- datExpr %>% t() 
# datMeta0 <- datMeta[c('AgeDeath', 'RIN1', 'PMI', 'Region', 'Sex')]
# datMeta0$Region <- as.numeric(datMeta0$Region)
# datMeta0$Sex <- as.numeric(datMeta0$Sex)
# lm.cov <- lm(dat ~ AgeDeath + RIN1 + PMI + Region + Sex, data=datMeta0) ##for sex specific, not including Race/Area
# dat.ori <- dat
# dat <- dat - as.matrix(datMeta0[c('AgeDeath', 'RIN1', 'PMI', 'Region', 'Sex')]) %*% lm.cov$coefficients[2:6,]

##
# enableWGCNAThreads()
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 3)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft$fitIndices[,5])

## module construction
# R2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
# power <- min(which(R2>.8))
power <- 4
net <- blockwiseModules(datExpr = dat, power = power,
                        TOMType = "signed", minModuleSize = 20,
                        reassignThreshold = 0, mergeCutHeight = 0.1,
                        numericLabels = TRUE, pamRespectsDendro = F,
                        saveTOMs = F, 
                        saveTOMFileBase = "catTOM",
                        verbose = 3)
# pdf('result_0703//wgcna_dendrogram_amp_ca.pdf', width=12)
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()


## module-trait cor
# datMeta1 <- data.frame(mdd=as.numeric(datMeta$PrimaryDx=="MDD"),
#                        ptsd=as.numeric(datMeta$PrimaryDx=="PTSD"))
datMeta1 <- data.frame(ptsd=as.numeric(datMeta$PrimaryDx=="PTSD"))
# datMeta1 <- data.frame(mdd=as.numeric(datMeta$PrimaryDx=="MDD"))
# datMeta1$mdd[datMeta$PrimaryDx=="PTSD"] <- NA
# datMeta1$ptsd[datMeta$PrimaryDx=="MDD"] <- NA
datMeta1 <- cbind(datMeta1, datMeta0)
moduleColors <- labels2colors(net$colors)
nGenes = ncol(dat);
MEs0 = moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# nSamples <- sum(!is.na(datMeta1$ptsd))
nSamples <- dim(datMeta1)[1]
moduleTraitCor = cor(MEs, datMeta1, use = "p");# Pearson correlation
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# pdf('result_0703/wgcna_traitcor_dl_cp.pdf')
# pdf('result_0703/wgcna_traitcor_dl_cm.pdf')
pdf('result_0703/wgcna_traitcor_sg_cp.pdf')
par(mar = c(5, 10, 3, 3));
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

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(dat, mydataP, use = "p"));
geneTraitSignificance = as.data.frame(cor(dat, datMeta1, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(datMeta1), sep="");
names(GSPvalue) = paste("p.GS.", names(datMeta1), sep="");
all <- data.frame(Geneid=datProbes$GN,
                  Protein=datProbes$EN,
                  module=moduleColors)
all <- cbind(all, trait=geneTraitSignificance)
mm <- sapply(1:dim(all)[1], function(x){geneModuleMembership[x, match(moduleColors[x], modNames)]})
all$DxmoduleMembership <- mm

##save
# filesave = "result_0703//WGCNA_CP_dl.RData"
# filesave =  "result_0703/WGCNA_CM_dl.RData"
# filesave = "result_0703/WGCNA_CM_sg.RData"
# save(dat, datProbes, datMeta, all, moduleTraitCor, moduleTraitPvalue, MEs, textMatrix, 
#      file=filesave)


### CSEA ====
# load(filesave)
load('../data/datProbes.RData')
names(all)[1] <- "Genename"
all$Geneid <- probes$ensembl_gene_id[match(all$Genename,probes$external_gene_name)]

load('../data/zhang.pSIout.RData')
colors = as.vector(all$module)
cell.p.zhang = matrix(NA, length(unique(colors)), 5);  rownames(cell.p.zhang) = unique(colors)
colnames(cell.p.zhang) = colnames(pSI.output)
for(mod in unique(colors)) {
  f = fisher.iteration(pSI.output, all$Geneid[all$module==mod], p.adjust = F)
  cell.p.zhang[mod,] = f$`0.05 - nominal`
}

# modSig <- rownames(moduleTraitPvalue)[moduleTraitPvalue[,1]<.05]
# modSig <- gsub('ME','',modSig)
# modSig <- modNames
modSig <- gsub('ME','',rownames(moduleTraitPvalue))
names(MEs) <- gsub('ME', '', names(MEs))
cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
# indSig <- match(modSig, rownames(cell.p.zhang.fdr))
# to_plot = cell.p.zhang.fdr[indSig,]
to_plot = cell.p.zhang.fdr
dendro.col = dendro.col.zhang
dendro.row= as.dendrogram(hclust(as.dist(1-bicor(MEs[modSig])), method="average"))

library(gplots)
library(WGCNA)
pdf('result_0703/wgcna_csea_sg_cp.pdf')
heatmap.2(-log10(to_plot+10^(-100)),
          col=blueWhiteRed(1000,1)[500:1000],
          scale="none",
          trace="none",
          cexRow = 0.8,
          cexCol = .8, 
          density.info = "none",
          colsep=0:7,
          rowsep=0:nrow(to_plot),
          sepcolor="grey",
          sepwidth=c(0.02,0.02),
          srtCol=45,
          offsetRow=0,
          offsetCol=-0.5,
          Rowv=dendro.row,
          Colv=dendro.col,
          key=T,
          key.xlab="-log10(P)",
          cellnote=signif(to_plot,1),
          notecex=.8,
          notecol="black",
          main='Enrichment')
dev.off()

##save
# filesave = "result_0703/WGCNA_DIA_dl_cp.RData"
# filesave = 'result_0703/WGCNA_CM_dl.RData'
filesave = 'result_0703/WGCNA_CP_sg.RData'
save(dat, datProbes, datMeta, all, to_plot, dendro.row, dendro.col, moduleTraitCor, moduleTraitPvalue, MEs, textMatrix, 
     file=filesave)


### write out genes for GO ====
## dl: cyan, yellow, salmon
load('result_0518/WGCNA_DIA_dl.RData')
genes <- all$Geneid[all$module=="salmon"]
write(genes %>% as.character, file = "result_0518/mod_salmon_dl.txt")

### WGCNA on combined dlPFC and sgPFC, CP/CM ====
library(magrittr)
library(stats)
library(WGCNA)
library(flashClust)
library(biomaRt)
library(corrplot)
library(ggplot2)
library(cqn)
library(sva)
library(limma)
library(statmod)
library(dplyr)
library(pSI)
library(pSI.data)
library(gplots)
options(stringsAsFactors = FALSE, digits = 3)

rm(list=ls())
load('DIA_new.RData')
dat <- datExpr %>% t()
gene.sel2 <- colSums(is.na(dat))==0
dat <- dat[,gene.sel2]
sam.sel <- datMeta$PrimaryDx!="PTSD" ## "MDD
demos <- datMeta[sam.sel,]
dat <- dat[sam.sel,]
suff <- datProbes[gene.sel2,]

datMeta <- demos[c("Sex","AgeDeath", "PMI", "Race", "Smoking", "RIN","Region")]
# datMeta$Dx.PTSD <- 0; datMeta$Dx.PTSD[demos$PrimaryDx=="PTSD"] <- 1
datMeta$Dx.MDD <- 0; datMeta$Dx.MDD[demos$PrimaryDx=="MDD"] <- 1
datMeta$Race <- 0; datMeta$Race[demos$Race=="W"] <- 1; datMeta$Race[demos$Race=="B"] <- -1 
datMeta$Smoking <- as.character(datMeta$Smoking); datMeta$Smoking[datMeta$Smoking=="Y"] <- 1; datMeta$Smoking[datMeta$Smoking=="N"] <- -1; datMeta$Smoking[datMeta$Smoking=="U"] <- 0
datMeta$Smoking <- as.numeric(datMeta$Smoking)
datMeta$Sex <- as.numeric(demos$Sex)
datMeta$Region <- as.numeric(datMeta$Region)
traitColors <- numbers2colors(datMeta, signed = FALSE)

# mod = model.matrix(~Dx.PTSD, data=datMeta)
mod = model.matrix(~Dx.MDD, data=datMeta)
datExpr.combat = ComBat(dat = as.matrix(t(dat)), batch = factor(datMeta$Region), mod = mod)
datExpr.preCombat = dat
dat = t(datExpr.combat)

##** regression model: reg2**
lm.cov <- lm(dat ~ Sex + AgeDeath + RIN + PMI + Race, data=datMeta) ##
dat.ori <- dat
dat <- dat - as.matrix(datMeta[c('Sex','AgeDeath', 'RIN', 'PMI', 'Race')]) %*% lm.cov$coefficients[2:6,]

##
disableWGCNAThreads()
powers <- 1:10
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 3)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft$fitIndices[,5])


## module construction
power <- 6
net <- blockwiseModules(dat, power = power,
                        TOMType = "signed", minModuleSize = 20, 
                        reassignThreshold = 0, mergeCutHeight = 0.1,
                        numericLabels = TRUE, pamRespectsDendro = F,
                        saveTOMs = F, saveTOMFileBase = "catTOM",
                        verbose = 3)

plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


## module-trait cor
moduleColors <- labels2colors(net$colors)
nGenes = ncol(dat);
MEs0 = moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
nSamples <- dim(datMeta)[1]
# moduleTraitCor = cor(MEs, datMeta, use = "p");# Pearson correlation
# datMeta1 <- datMeta[c("Race","Smoking","Dx.MDD","Dx.PTSD","PrimaryDx","BA11","BA24","BA25","BA9")]
# datMeta1 <- datMeta[c("Race","Smoking","Dx.MDD","Dx.PTSD","PrimaryDx")]
# datMeta1 <- datMeta[c('Dx.PTSD','Smoking')]
datMeta1 <- datMeta[c('Dx.MDD','Smoking')]
moduleTraitCor = cor(MEs, datMeta1, use = "p") #
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datMeta1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(dat, mydataP, use = "p"));
geneTraitSignificance = as.data.frame(cor(dat, datMeta, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(datMeta), sep="");
names(GSPvalue) = paste("p.GS.", names(datMeta), sep="");
all <- data.frame(PN=datProbes$PN,
                  Genename=datProbes$GN,
                  module=moduleColors)
all <- cbind(all, traitSignificance=geneTraitSignificance)
mm <- sapply(1:dim(all)[1], function(x){geneModuleMembership[x, match(moduleColors[x], modNames)]})
all$DxmoduleMembership <- mm

save(all, moduleTraitCor, moduleTraitPvalue, dat, datMeta, file="result_0703/WGCNA_CM_dlsg.RData")
pdf('result_0703/wgcna_traitcor_cm_dlsg.pdf')
par(mar = c(5, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datMeta1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


### consistency between CP/CM modules: enrichment ====
library(ggplot2)
library(dplyr)
rm(list=ls())
load('result_0703/WGCNA_CP_dlsg.RData')
all.p <- all
load('result_0703/WGCNA_CM_dlsg.RData')
all.t <- all

all <- merge(all.p[2:3], all.t[2:3], by="Genename")
table(all[c("module.x","module.y")]) %>% View #%>%heatmap
df <- table(all[c("module.x","module.y")]) %>% as.data.frame
df$ES <- df$Freq*dim(all)[1]/as.numeric(table(all$module.x)[df$module.x]*table(all$module.y)[df$module.y])
names(df)[1:2] <- c("module.protein","module.RNA")
max.mp <- sapply(levels(df$module.protein), function(m) max(df$ES[df$module.protein==m])) %>% sort(decreasing = T)
max.mt <- sapply(levels(df$module.RNA), function(m) max(df$ES[df$module.RNA==m])) %>% sort
df$module.protein <- factor(df$module.protein, levels=names(max.mp))
df$module.RNA <- factor(df$module.RNA, levels=names(max.mt))

ggplot(df, aes(module.protein, module.RNA, fill=log1p(Freq))) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_classic()

pdf('result_0703/wgcna_module_cp_cm_dlsg.pdf', width=8, height=6)
ggplot(df[!(df$module.protein=="grey"|df$module.RNA=="grey"),], 
       aes(module.protein, module.RNA, fill=ES)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  xlab('module.PTSD') + ylab('module.MDD') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### comparing RNA/protein expression: select the same individuals ====
rm(list=ls())
load('result_upmc/dat_upmc_qnorm_CMP.RData')
datExpr.t <- datExpr; datMeta.t <- datMeta; datProbes.t <- datProbes
load('DIA_new.RData')

datExpr <- as.data.frame(datExpr)
sel_sam <- which(datMeta.t$Region==9)
datExpr.t <- datExpr.t[,sel_sam]
datMeta.t <- datMeta.t[sel_sam,]
datExpr.t$Geneid <- datProbes.t$Genename
datExpr$Geneid <- datProbes$GN
dat <- left_join(datExpr[,c(which(datMeta$Region=="dlPFC"),ncol(datExpr))], datExpr.t, by="Geneid")
mean_prot <- apply(dat[,1:(grep("Geneid",names(dat))-1)], 1, mean)
mean_rna <- apply(dat[,(grep("Geneid",names(dat))+1):ncol(dat)], 1, mean)
lmm <- lm(mean_rna ~ mean_prot) %>% summary %>% .[['coefficients']] %>% .[,1]
plot(mean_prot, mean_rna, pch=16, col=rgb(0,0,0,.3))
sel.top <- match(genes.top4, dat$Geneid)
sel.bot <- match(genes.bot2, dat$Geneid)
points(mean_prot[sel.top], mean_rna[sel.top], pch=16, col="green")
points(mean_prot[sel.bot], mean_rna[sel.bot], pch=16, col="red")
abline(lmm)

## individual wise
brnums <- intersect(datMeta$BrNum, datMeta.t$BrNum)
r2 <- sapply(brnums, function(b) lm(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])] ~
                                      dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)]) %>% summary %>% .[['adj.r.squared']])
hist(r2, breaks=20)
b = brnums[1]
res1 <- lm(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])] ~ dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)])
res2 <- MASS::rlm(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])] ~ dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)])
plot(dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)], 
     dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])],
     col=rgb(0,0,0,0.2), pch=16)
abline(res1$coefficients)
## only a very slight change in slope
r2 <- sapply(brnums, function(b) lm(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])] ~
                                      dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)]) %>% summary %>% .[['coefficients']] %>% .[,1])
r2rb <- sapply(brnums, function(b) MASS::rlm(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])] ~
                                               dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)]) %>% summary %>% .[['coefficients']] %>% .[,1])
hist(r2[2,], breaks=20)
hist(r2rb[2,], breaks=20)
corsp <- sapply(brnums, function(b) cor(dat[,match(b,datMeta$BrNum[datMeta$Region=="dlPFC"])],
                                        dat[,grep("Geneid",names(dat))+match(b, datMeta.t$BrNum)], 
                                        use="pairwise.complete.obs", method="pearson"))
hist(corsp, breaks=20, xlab="Spearman corr")

## gene wise
dat1 <- dat[,1:(grep("Geneid",names(dat))-1)] ## protein
dat2 <- dat[,(grep("Geneid",names(dat))+1):ncol(dat)] ## RNA
dat2 <- dat2[,match(datMeta$BrNum[datMeta$Region=="dlPFC"],datMeta.t$BrNum)]
gene.sel <- !is.na(rowSums(dat2)) & !is.na(rowSums(dat1))
dat1 <- dat1[gene.sel,]
dat2 <- dat2[gene.sel,]
corsp <- cor(dat1, dat2, use="pairwise.complete.obs", method="pearson") ## individual wise
hist(diag(corsp), breaks=20, xlab="Pearson corr")
heatmap(corsp, symm=T, scale = "none")
corsp <- cor(t(dat1), t(dat2), use="pairwise.complete.obs", method="pearson") ## gene wise
hist(diag(corsp), breaks=20, xlab="Pearson corr")
set.seed(1116)
subsam <- sample(1:nrow(corsp), 500)
heatmap(corsp[subsam,subsam], symm=T, scale="none")


### GO of top correlated genes
genes <- dat$Geneid[gene.sel]
genes.top4 <- genes[diag(corsp)>.4]
genes.bot2 <- genes[diag(corsp)< -0.2]
library(gprofiler2)
reportGo <- function(query){
  go <- gost(query, organism = "hsapiens", significant = T, evcodes=T)$result
  go <- go[order(go$p_value, decreasing = F),]
  go[,c(1:6,9:11,16)]
}
# gos <- reportGo(genes.top4)
gos <- reportGo(genes.bot2)
View(gos[c(3,8,9)])


### are they homogenous, RNA and protein expression
dat1n <- apply(dat1, 1, scale) %>% t()
dat2n <- apply(dat2, 1, scale) %>% t()
dat3n <- cbind(dat1n,dat2n)
dat31 <- dat3n[sample(1:nrow(dat3n),100),]
heatmap(cor(t(dat31) %>% as.matrix))


### comparing RNA/protein expression: corresponding UPMC samples only ====
rm(list=ls())
load('proteomics/DIA_new.RData')
datExpr.p <- datExpr; datMeta.p <- datMeta; datProbes.p <- datProbes
load('data/dat13_ulval_keep3.RData')
sel <- meta13uv_keep$BrNum %in% (datMeta$BrNum[datMeta$Region=="dlPFC"] ) 
datExpr.t <- dat13uv_keep[,6+which(sel)]
datMeta.t <- meta13uv_keep[sel,]
fpkm.t <- fpkm13uv_keep[,6+which(sel)]
save(datExpr.t, datMeta.t, suff_keep, fpkm.t, file="data/RNA_upmc_prot_dl.RData")

### DESeq2 with the 52 UPMC samples ====
rm(list=ls())
load('RNA_upmc_prot_dl.RData')

library(DESeq2)
cntdat <- datExpr.t
demos <- datMeta.t
suff <- suff_keep
demos$RIN2 <- demos$RIN1^2
demos$Area <- demos$Region
demos$AgeDeath <- as.numeric(demos$AgeDeath)
demos$PMI <- as.numeric(demos$PMI)
demos$Race <- factor(demos$Race)
rownames(demos) <- NULL

bysex = T
if(bysex){
  indices <- list(which(demos$Area==11 & demos$Sex=="F"),
                  which(demos$Area==25 & demos$Sex=="F"),
                  which(demos$Area==24 & demos$Sex=="F"),
                  which(demos$Area==9  & demos$Sex=="F"),
                  which(demos$Area==11 & demos$Sex=="M"),
                  which(demos$Area==25 & demos$Sex=="M"),
                  which(demos$Area==24 & demos$Sex=="M"),
                  which(demos$Area==9  & demos$Sex=="M")) 
  files <- paste(rep(c("F","M"),each=4), c(11,25,24,9), sep="")
} else{
  indices <- list(which(demos$Area==11),
                  which(demos$Area==25),
                  which(demos$Area==24),
                  which(demos$Area==9))
  files <- paste0("A", c(11,25,24,9))
}
alpha <- .5 ##by default .5
for (i in 1:length(indices)){
  index <- indices[[i]] 
  genes <- (rowSums(cntdat[,index]) > dim(cntdat)[2]*alpha)
  seldat <- cntdat[genes,index]
  if(bysex){
    dds <- DESeqDataSetFromMatrix(countData = seldat,
                                  colData = demos[index,],
                                  # design= ~ PrimaryDx + AgeDeath + RIN1) ##model 3
                                  design = ~ PrimaryDx + AgeDeath + RIN1 + PMI + Race) ##model 4
  }else   dds <- DESeqDataSetFromMatrix(countData = seldat,
                                        colData = demos[index,],
                                        design= ~ PrimaryDx + AgeDeath + RIN1 + Sex + PMI + Race)
  
  dds <- DESeq(dds)
  for (j in 1:3){
    if (j == 1) {
      res <- results(dds, contrast=c("PrimaryDx", "MDD", "Control"))
    } else if(j ==2) {
      res <- results(dds, contrast=c("PrimaryDx", "PTSD", "Control"))
    } else res <- results(dds, contrast=c("PrimaryDx", "PTSD", "MDD"))
    
    dd <- res@listData %>% as.data.frame
    dd <- cbind(dd, data.frame(Geneid=suff$Geneid[genes], Genename=suff$Genename[genes]))
    dsort <- dd[order(dd$padj),]
    # dsort$padj <- p.adjust(dsort$pvalue, method = "fdr") ##somehow padjust failed for F25 calculation
    dsort <- dsort[!is.na(dsort$padj),]
    # dsort <- dsort[!is.na(dsort$pvalue),]
    if (j == 1) {
      sig01 <- dsort
    } else if(j == 2){
      sig02 <- dsort
    } else sig12 <- dsort
  }
  fname <- paste0("result_upmc/deseq_", files[i], "_m4_keep3.RData")
  save(sig01, sig02, sig12, file=fname)
}

### utmost results for multiple traits ####
rm(list=ls())
## read all utmost
ad <- read.table('result_gwas/merged_test_results_MDD.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; ad_sig <- ad[1:100,] ## only the top 100
# ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
# sum(ad$p_adj<.05)
# ad_sig <- ad[ad$p_adj<.05,]
# ad_sig <- ad_sig[order(ad_sig$p_value),]
mdd <- ad_sig
ad <- read.table('result_gwas/merged_test_results_SCZ.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; ad_sig <- ad[1:100,] ## only the top 100
# ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
# sum(ad$p_adj<.05)
# ad_sig <- ad[ad$p_adj<.05,]
# ad_sig <- ad_sig[order(ad_sig$p_value),]
scz <- ad_sig
ad <- read.table('result_gwas/merged_test_results_BD.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; ad_sig <- ad[1:100,] ## only the top 100
# ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
# sum(ad$p_adj<.05)
# ad_sig <- ad[ad$p_adj<.05,]
# ad_sig <- ad_sig[order(ad_sig$p_value),]
bpd <- ad_sig
ad <- read.table('../results/gwas/merged_test_results_all.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; ad_sig <- ad[1:100,] ## only the top 100
# ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
# sum(ad$p_adj<.05)
# ad_sig <- ad[ad$p_adj<.05,]
# ad_sig <- ad_sig[order(ad_sig$p_value),]
psd <- ad_sig

## put them together
azd <- read.table('result_gwas/merged_test_results_AD.txt', header=T)
azd <- azd[!is.na(azd$test_score),]
azd <- azd[order(azd$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
azd$p_adj <- p.adjust(azd$p_value, method="bonferroni")
ad <- read.table('result_gwas/merged_test_results_MDD.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
mdd <- ad
ad <- read.table('result_gwas/merged_test_results_SCZ.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
scz <- ad
ad <- read.table('result_gwas/merged_test_results_BD.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
bpd <- ad
ad <- read.table('../results/gwas/merged_test_results_all.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
psd <- ad
ad <- read.table('result_gwas/merged_test_results_ASD.txt', header=T)
ad <- ad[!is.na(ad$test_score),]
ad <- ad[order(ad$p_value),]; #ad_sig <- ad[1:100,] ## only the top 100
ad$p_adj <- p.adjust(ad$p_value, method="bonferroni")
asd <- ad
save(asd, azd, mdd, psd, scz, bpd, file="utmost_6traits.RData")

####### END ##############
