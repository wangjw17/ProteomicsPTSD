#######################################################
######### miRNA analysis of proteomics project ########
#######################################################

### STAR params for miRNAseq mapping ####
## https://github.com/sorenar/miRNA-seq-adapters/blob/master/cutadapt_set1.sh
## https://www.encodeproject.org/pipelines/ENCPL444CYA/ 
## recommended:
--runThreadN 16
--alignEndsType EndToEnd ## default: local
--alignSJDBoverhangMin 1000 ## default: 3
--alignIntronMax 1 ## default: 0
--quantMode TranscriptomeSAM GeneCounts ## default: -
--outSAMtype BAM SortedByCoordinate ## default: SAM Unsorted/SortedByCoordinate
--outSAMunmapped Within ## default: none
--outReadsUnmapped Fastx ## default: none
--outFilterScoreMinOverLread 0 ## default: 0.66
--outFilterMismatchNmax 1 ## default 10
--outFilterMultimapNmax 10 ## default
--outFilterMultimapScoreRange 0 ## default: 1
--outFilterMatchNmin 16 ## default: 0
--outFilterMatchNminOverLread 0 ## default: 0.66
--outWigType wiggle ## default: none
--outWigStrand Stranded ## default
--outWigNorm RPM ## default
## modified
params=' --runThreadN 12
--sjdbGTFfile /home/jw2372/scratch60/ptsd/hsp_genome/Homo_sapiens.GRCh38.104.gtf
--alignSJDBoverhangMin 1000
--outSAMtype BAM SortedByCoordinate
--outFilterMultimapScoreRange 0
--outFilterMatchNmin 16
--outFilterMatchNminOverLread 0.3
--outFilterScoreMinOverLread 0.3 
'
genomeDir=/home/jw2372/scratch60/ptsd/hsp_genome
STAR --genomeDir $genomeDir --readFilesIn RTMIR040_033_CAGGCG_L001_R1_001.fastq $params


### read in miRNA-seq reads ####
library(stringr)
rm(list=ls())
# dat <- read.table('rawdata/newcounts_all1.5.txt', header=T) ## mapping trial1.4; counting 1.5; -t exon -g transcript_id -s 1 -O

datProbes <- dat[1:6]
datProbes$chr <- datProbes$Chr %>% as.character %>% strsplit(split=";") %>% sapply("[",1)
datExpr <- dat[-(1:6)]
rm(dat)
samples <- names(datExpr) %>% strsplit(split="\\.") %>% sapply("[",1) %>% strsplit(split="_") %>% sapply("[",2)
samples0 <- names(datExpr) %>% strsplit(split="\\.") %>% sapply("[",1)
indnum <- names(datExpr) %>% strsplit(split="\\.") %>% sapply("[",1) %>% strsplit(split="_") %>% sapply("[",3) %>% as.numeric
names(datExpr) <- samples
datMeta <- data.frame(Sample0=samples0,
                      Sample=samples,
                      IndNum=indnum)
# BrNum=str_extract(samples, pattern="[0-9]+[U,C,M,P]") %>% gsub("U|C|M|P","",.),
# PrimaryDx=str_extract(samples, pattern="[C,M,P]"),
# Region=str_extract(samples, pattern="[C,M,P][0-9]+") %>% gsub("C|M|P","",.))

# log <- read.table('rawdata/log_summary_all1.5.txt', header=T) %>% t()

colnames(log) <- log[1,] %>% as.character
log <- log[-1,] %>% as.data.frame
log$`Number of input reads` <- log$`Number of input reads` %>% as.character %>% as.numeric
log$`Uniquely mapped reads %` <- log$`Uniquely mapped reads %` %>% as.character %>% gsub("%","",.) %>% as.numeric
log$Assigned <- log$Assigned %>% as.character %>% as.numeric
rownames(log) <- datMeta$Sample

log$`Number of input reads` %>% log10 %>% hist(breaks=30,main="Histogram of total reads", xlab="log10(counts)")
log$`Uniquely mapped reads %` %>% hist(main="Histrogram of Uniquely mapped reads %")
log$Assigned %>% log10 %>% hist(breaks=30,main="Histogram of Assigned reads log10")
(log$Assigned/log$`Number of input reads`) %>% hist(main="Histogram of Turn out rate")
log$`% of reads unmapped: too short` %>% gsub("%","",.) %>% as.numeric %>% hist(main="Histrogram of Unmapped too short %")


gene0 <- rowSums(datExpr)>0 & datProbes$chr %in% c(1:22, "X","Y","MT")
datExpr <- datExpr[gene0,]; datProbes <- datProbes[gene0,]
datProbes$chr <- factor(datProbes$chr)
dim(datExpr)
# save(datExpr, datMeta, datProbes, log, file="rawdata/dat_upmc_mirna_map1.5.RData")


# ref <- xlsx::read.xlsx('rawdata/Copy of RQ14148-summary.xlsx', sheetIndex=1)
# ref$BA <- rep(c(25,9), each=66)
# ref$label.dx <- str_extract(ref$NA., pattern="[C,M,P]")
# ref$Br.[ref$Br.==1988] <- 1488
# ref$Br.[ref$Br.==1912] <- 1412
# ref$Br.[ref$Br.==1204] <- 1201
# ref$label <- paste0(ref$Br.,"U",ref$label.dx,ref$BA)
# save(ref, file="rawdata/reference_mirna.RData")
load('rawdata/reference_mirna.RData')
datMeta$BrNum <- ref$Br.[match(datMeta$Sample,ref$Sample.ID)]
datMeta$Dx <- ref$NA.[match(datMeta$Sample,ref$Sample.ID)]
datMeta$Sample.bulk <- ref$label[match(datMeta$Sample,ref$Sample.ID)]


meta <- read.csv('../data/Sample_metadata.csv')
ref$Br. %in% meta$BrNum
datMeta0 <- datMeta
datMeta <- meta[match(datMeta$Sample.bulk, meta$Sample %>% gsub("r","",.)),]
datMeta$Sample.mir <- datMeta0$Sample
# names(datExpr) <- datMeta$Sample
## some of the samples are not correct: 008,033 not present
# sum(is.na(datMeta$BrNum))
# ref[ref$Sample.ID %in% paste0("RTMIR0",c("08","33")),] ##1204UC25, 1912UP25



## for transcript level
datProbes0 <- datProbes[c(1,6,7)]
load('../data/datProbes_transcript.RData')
names(datProbes0)[1] <- "TranscriptID"
datProbes0$Geneid <- probes$Gene.stable.ID[match(datProbes0$TranscriptID, probes$Transcript.stable.ID)]
datProbes0$Genename <- probes$Gene.name[match(datProbes0$TranscriptID, probes$Transcript.stable.ID)]
datProbes0$TranscriptLen <- probes$Transcript.length..including.UTRs.and.CDS.[match(datProbes0$TranscriptID, probes$Transcript.stable.ID)]
datProbes0$TranscriptType <- probes$Transcript.type[match(datProbes0$TranscriptID, probes$Transcript.stable.ID)]
datProbes0$TranscriptName <- probes$Transcript.name[match(datProbes0$TranscriptID, probes$Transcript.stable.ID)]
datProbes <- datProbes0

## for gene level
datProbes0 <- datProbes[c(1,6,7)]
load('../data/datProbes.RData')
datProbes0$Genename <- probes$external_gene_name[match(datProbes0$Geneid, probes$ensembl_gene_id)]
datProbes0$Biotype <- probes$gene_biotype[match(datProbes0$Geneid, probes$ensembl_gene_id)]
datProbes <- datProbes0

save(datExpr, datProbes, datMeta, log, file="rawdata/dat_upmc_mirna_t0m1.4.RData")

### DE miRNA analysis ####
rm(list=ls())
# load('rawdata/dat_upmc_mirna_trial4_etO.RData')
# load('rawdata/dat_upmc_mirna_map1.5.RData')
load('rawdata/dat_upmc_mirna_all1.5.RData')
datMeta$region <- factor(datMeta$region)
gene.mir <- datProbes$TranscriptType == "miRNA"
datExpr <- datExpr[gene.mir,]
datProbes <- datProbes[gene.mir,]
table(datMeta[c('PrimaryDx','Sex','region')])

library(DESeq2)
# sel.sam <- which(datMeta$region %in% c("dlPFC"))
sel.sam <- 1:nrow(datMeta)
sel.gene <- which(datProbes$TranscriptType=="miRNA")
cntdat <- datExpr[sel.gene,sel.sam]
demos <- datMeta[sel.sam,]
suff <- datProbes[sel.gene,]

demos$AgeDeath <- as.numeric(demos$AgeDeath)
demos$PMI <- as.numeric(demos$PMI)
demos$Race <- factor(demos$Race)
demos$AgeDeath <- as.numeric(demos$AgeDeath)
rownames(demos) <- NULL

bysex = F
if(bysex){
  indices <- list(which(demos$region=='dlPFC' & demos$Sex=="F"),
                  which(demos$region=='dlPFC' & demos$Sex=="M"))
  files <- paste(c("F","M"), 9, sep="")
} else{
  indices <- list(which(demos$region=='dlPFC'))
  files <- paste0("A", c(9))
}
alpha <- .5 ##by default .5
for (i in 1:length(indices)){
  index <- indices[[i]] 
  genes <- (rowSums(cntdat[,index]) > dim(cntdat)[2]*alpha)
  seldat <- cntdat[genes,index]
  dim(seldat)
  if(bysex){
    dds <- DESeqDataSetFromMatrix(countData = seldat,
                                  colData = demos[index,],
                                  design= ~ PrimaryDx + AgeDeath + RIN) ##model 3
    # design = ~ PrimaryDx + AgeDeath + RIN + PMI + Race) ##model 4
  }else{
    dds <- DESeqDataSetFromMatrix(countData = seldat,
                                  colData = demos[index,],
                                  design= ~ PrimaryDx + AgeDeath + RIN + Sex) ##model 3
    # design= ~ PrimaryDx + AgeDeath + RIN + Sex + PMI + Race) ##model 4
    
  }
  
  dds <- DESeq(dds)
  for (j in 1:3){
    if (j == 1) {
      res <- results(dds, contrast=c("PrimaryDx", "MDD", "Control"))
    } else if(j ==2) {
      res <- results(dds, contrast=c("PrimaryDx", "PTSD", "Control"))
    } else res <- results(dds, contrast=c("PrimaryDx", "PTSD", "MDD"))
    
    dd <- res@listData %>% as.data.frame
    dd <- cbind(dd, data.frame(Geneid=suff$Geneid[genes],
                               Chr=suff$chr[genes],
                               Transcriptid=suff$Genename[genes]))
    dsort <- dd[order(dd$pvalue),]
    # dsort$padj <- p.adjust(dsort$pvalue, method = "fdr") ##somehow padjust failed for F25 calculation
    dsort <- dsort[!is.na(dsort$pvalue),]
    if (j == 1) {
      sig01 <- dsort
    } else if(j == 2){
      sig02 <- dsort
    } else sig12 <- dsort
  }
  # fname <- paste0("result_mirna/deseq_", files[i], "_m3_ulval_transcript_lncproc.RData")
  fname <- paste0("result_mirna/deseq_", files[i], "_m3_upmc_mirna_all1.5.RData")
  save(sig01, sig02, sig12, file=fname)
}

## summary
rm(list=ls())
# labs <- paste0(rep(c('A','F','M'),each=4),c(9,11,24,25))
labs <- paste0(rep(c('A','F','M'),each=2),c(9,25))
for (lab in labs){
  # file <- paste0('result_upmc/deseq_',lab,'_m4_re.RData')
  # file <- paste0('result_mirna/deseq_',lab,'_m3_upmc_mirna_etO.RData')
  # file <- paste0('result_mirna/deseq_',lab,'_m3_upmc_mirna_map1.5.RData')
  file <- paste0('result_mirna/deseq_',lab,'_m3_upmc_mirna_all1.5.RData')
  if (file.exists(file)){
    load(file)
    lab %>% print
    dim(sig01) %>% print
    sum(sig01$pvalue<.05) %>% print
    sum(sig01$padj<.05,na.rm=T) %>% print
    sum(sig02$pvalue<.05) %>% print
    sum(sig02$padj<.05,na.rm=T) %>% print
  }
}

load('result_mirna/deseq_A9_m3_upmc_mirna_etO.RData')
load('rawdata/hsa_probes.RData')
sig02$Name <- hsa_gtf$gene_name[match(sig02$Transcriptid %>% as.character,hsa_gtf$transcript_id)]
sig02$Name2 <- hsa_gtf$gene_name[match(sig02$Geneid,hsa_gtf$gene_id)]
View(sig02)

### comparing between dlPFC and sgPFC, MDD and PTSD, FvsM ####
rm(list=ls())
# load('result_mirna/deseq_A9_m3_upmc_mirna_map1.5.RData')
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
sig01.d <- sig01; sig02.d <- sig02
load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
sig01.s <- sig01; sig02.s <- sig02

## comparing brain regions
sig01 <- merge(sig01.d, sig01.s, by="Geneid")
sig01 <- sig01[sig01$Chr.x!="X" & !is.na(sig01$padj.x) & !is.na(sig01$padj.y),]
lm.01 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig01) %>% summary
plot(sig01$log2FoldChange.x, sig01$log2FoldChange.y, pch=16, xlab="dlPFC", ylab="sgPFC", 
     main=paste0("MDD, R2=",lm.01 %>% .$adj.r.squared %>% round(digits=3)))
sel.sig <- apply(sig01[c('padj.x','padj.y')],1,min) < 0.05
points(sig01$log2FoldChange.x[sel.sig], sig01$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.01$coefficients[,1], col="red", lty=2)
sig02 <- merge(sig02.d, sig02.s, by="Geneid")
sig02 <- sig02[sig02$Chr.x!="X" & !is.na(sig02$padj.x) & !is.na(sig02$padj.y),]
lm.02 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig02) %>% summary 
sel.sig <- apply(sig02[c('padj.x','padj.y')],1,min) < 0.05
plot(sig02$log2FoldChange.x, sig02$log2FoldChange.y, pch=16, xlab="dlPFC", ylab="sgPFC", 
     main=paste0("PTSD, R2=",lm.02 %>% .$adj.r.squared %>% round(digits=3)))
points(sig02$log2FoldChange.x[sel.sig], sig02$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.02$coefficients[,1], col="red", lty=2)

## comparing MDD and PTSD
sig01 <- merge(sig01.d, sig02.d, by="Geneid")
sig01 <- sig01[sig01$Chr.x!="X" & !is.na(sig01$padj.x) & !is.na(sig01$padj.y),]
lm.01 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig01) %>% summary
sel.sig <- apply(sig01[c('padj.x','padj.y')],1,min) < 0.05
plot(sig01$log2FoldChange.x, sig01$log2FoldChange.y, pch=16, xlab="MDD", ylab="PTSD", 
     main=paste0("dlPFC, R2=",lm.01 %>% .$adj.r.squared %>% round(digits=3)))
points(sig01$log2FoldChange.x[sel.sig], sig01$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.01$coefficients[,1], col="red", lty=2)
sig02 <- merge(sig01.s, sig02.s, by="Geneid")
sig02 <- sig02[sig02$Chr.x!="X" & !is.na(sig02$padj.x) & !is.na(sig02$padj.y),]
lm.02 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig02) %>% summary 
sel.sig <- apply(sig02[c('padj.x','padj.y')],1,min) < 0.05
plot(sig02$log2FoldChange.x, sig02$log2FoldChange.y, pch=16, xlab="MDD", ylab="PTSD", 
     main=paste0("sgPFC, R2=",lm.02 %>% .$adj.r.squared %>% round(digits=3)))
points(sig02$log2FoldChange.x[sel.sig], sig02$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.02$coefficients[,1], col="red", lty=2)

## comparing sex overlap
rm(list=ls())
# load('result_mirna/deseq_F25_m3_upmc_mirna_map1.5.RData')
load('result_mirna/deseq_F9_m3_upmc_mirna_all1.5.RData')
sig01.d <- sig01; sig02.d <- sig02
# load('result_mirna/deseq_M25_m3_upmc_mirna_map1.5.RData')
load('result_mirna/deseq_M9_m3_upmc_mirna_all1.5.RData')
sig01.s <- sig01; sig02.s <- sig02
sig01 <- merge(sig01.d, sig01.s, by="Geneid")
sig01 <- sig01[sig01$Chr.x!="X" & !is.na(sig01$padj.x) & !is.na(sig01$padj.y),]
lm.01 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig01) %>% summary
sel.sig <- apply(sig01[c('padj.x','padj.y')],1,min) < 0.05
plot(sig01$log2FoldChange.x, sig01$log2FoldChange.y, pch=16, xlab="Female", ylab="Male", 
     main=paste0("dlPFC MDD, R2=",lm.01 %>% .$adj.r.squared %>% round(digits=3)))
points(sig01$log2FoldChange.x[sel.sig], sig01$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.01$coefficients[,1], col="red", lty=2)

sig02 <- merge(sig02.d, sig02.s, by="Geneid")
sig02 <- sig02[sig02$Chr.x!="X" & !is.na(sig02$padj.x) & !is.na(sig02$padj.y),]
lm.02 <- lm(log2FoldChange.x~log2FoldChange.y, data=sig02) %>% summary 
sel.sig <- apply(sig02[c('padj.x','padj.y')],1,min) < 0.05
plot(sig02$log2FoldChange.x, sig02$log2FoldChange.y, pch=16, xlab="Female", ylab="Male", 
     main=paste0("dlPFC PTSD, R2=",lm.02 %>% .$adj.r.squared %>% round(digits=3)))
points(sig02$log2FoldChange.x[sel.sig], sig02$log2FoldChange.y[sel.sig], pch=16, col="red")
abline(0,1, lty=2)
abline(lm.02$coefficients[,1], col="red", lty=2)


### compare miRNA with proteins ####
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
# load('../data/ulval_re.RData') ## miRNA data from bulk RNAseq with Ensembl reference feature exon
# sel.sam <- datMeta$region=='dlPFC'
# sel.gene <- which(datProbes$gene_biotype=="miRNA")
# load('rawdata/ulval_re_mi_map2.RData') ## miRNA data from bulk RNAseq with MirBase reference
# load('rawdata/dat_upmc_mirna_map1.5.RData')
load('rawdata/dat_upmc_mirna_all1.5.RData')

# ## normalization on selected miRNAs (option 1)
# sel.sam <- datMeta$region=="dlPFC"
# sel.gene <- datProbes$TranscriptType=="miRNA"
# datExpr <- datExpr[sel.gene,sel.sam]; datProbes <- datProbes[sel.gene,]; datMeta <- datMeta[sel.sam,]
# datExpr.t <- datExpr; datProbes.t <- datProbes; datMeta.t <- datMeta
# sel.gene2 <- rowSums(datExpr.t) > ncol(datExpr.t) * 0.5
# datExpr.t <- datExpr.t[sel.gene2,]; datProbes.t <- datProbes.t[sel.gene2,]
# datFPKM.t <- getFKPM(datExpr.t, trim=0, suffix=F, lengths=datProbes.t$Length)
# rownames(datFPKM.t) <- datProbes.t$Geneid

## filter after normalization (option 2)
sel.sam <- datMeta$region=="dlPFC"
datExpr <- datExpr[,sel.sam]; datMeta <- datMeta[sel.sam,]
datFPKM.t <- getFKPM(datExpr, trim=0, suffix=F, lengths=datProbes$Length)
sel.gene <- datProbes$TranscriptType=="miRNA"
datFPKM.t <- datFPKM.t[sel.gene,]; datProbes.t <- datProbes[sel.gene,]; datMeta.t <- datMeta
datExpr.t <- datExpr[sel.gene,]
sel.gene2 <- rowSums(datExpr.t) > ncol(datExpr.t) * 0.5
datFPKM.t <- datFPKM.t[sel.gene2,]; datProbes.t <- datProbes.t[sel.gene2,]
rownames(datFPKM.t) <- datProbes.t$Geneid


## Dx specificity
sel <- datMeta.t$PrimaryDx != "MDD"
datFPKM.t <- datFPKM.t[,sel]; datMeta.t <- datMeta.t[sel,]

### proteins
load('DIA_new.RData')
datExpr.p <- datExpr; datMeta.p <- datMeta; datProbes.p <- datProbes
sel.sam <- datMeta.p$Region=='dlPFC'
datExpr.p <- datExpr.p[,sel.sam]; datMeta.p <- datMeta.p[sel.sam,]
datFPKM.t1 <- datFPKM.t[,match(datMeta.p$BrNum,datMeta.t$BrNum)]
rownames(datExpr.p) <- datProbes.p$EN_short
cors <- cor(t(datFPKM.t1),t(datExpr.p), use = "complete.obs")
cors %>% abs %>% apply(1,max) %>% hist(main="Max correlation for each miRNA", breaks=50)
cors %>% abs %>% apply(1,median) %>% hist(main="Median correlation for each miRNA", breaks=50)
df <- melt(cors)
df$pvalue <- sapply(colnames(cors), function(g2){ 
  sapply(1:nrow(cors), function(g1){cor.test(datFPKM.t1[g1,],datExpr.p[g2,])$p.value})}) %>% melt %>% .$value
# df$padj <- p.adjust(df$pvalue, method="fdr")
# sum(df$padj<.05)
# load('rawdata/hsa_probes.RData')
# df$Name <- hsa_mirbase$Name[match(df$Var1,hsa_mirbase$Alias)]
# df[df$Var1=="MIR128-1",] %>% .[order(abs(.$value), decreasing = T),] %>% View
# df[df$Var2=="VIAAT",] %>% .[order(abs(.$value), decreasing = T),] %>% View
# df[order(df$pvalue, decreasing=F),][1:100,] %>% View
# df[df$padj<.05,] %>% View
# df$Var2[df$padj<.05] %>% factor %>% table %>% sort(decreasing=T)
# df_sig <- df[df$pvalue<.05,]
# save(df_sig, file="result_mirna/cor_prot_dlpfc_sig.RData")

## combine with protein results
df$miRNA <- datProbes.t$Genename[match(df$Var1,datProbes.t$Geneid)]
load('result_mirna/deseq_A9_m3_upmc_mirna_all1.5.RData')
# load('result_mirna/deseq_A25_m3_upmc_mirna_map1.5.RData')
df$PTSD.padj <- sig02$padj[match(df$Var1, sig02$Geneid)]
df$PTSD.pvalue <- sig02$pvalue[match(df$Var1, sig02$Geneid)]
df$MDD.padj <- sig01$padj[match(df$Var1, sig01$Geneid)]
df$MDD.pvalue <- sig01$pvalue[match(df$Var1, sig01$Geneid)]
# df[which(df$padj<0.05),c(2:8)] %>% View
# df[which(df$MDD.pvalue<0.05 & df$padj<0.05),c(2:6,9:10)] %>% View
# df[which(df$padj<.05),c(2:7)] %>% View
# subset(df, PTSD.pvalue < .01 & pvalue < .01) %>% dim
# subset(df, MDD.pvalue < .01 & pvalue < .01) %>% dim
# subset(df, miRNA=="MIR218-1" & pvalue < .01 & PTSD.pvalue < .01) %>% View


### summarize MDD and PTSD DEGs for each region ####
degs <- list()
for (area in c("9","11","24","25")){
  load(paste0('dat_ulval/deseq_A', area, '_m3_ulval_keep3.RData'))
  sig01$Geneid <- as.character(sig01$Geneid); sig02$Geneid <- as.character(sig02$Geneid)
  sig01$Genename <- as.character(sig01$Genename); sig02$Genename <- as.character(sig02$Genename)
  sig01$Genename[sig01$Genename==""] <- sig01$Geneid[sig01$Genename==""]
  sig02$Genename[sig02$Genename==""] <- sig02$Geneid[sig02$Genename==""]
  ## not ranekd
  # degs[[area]] <- unique(intersect(sig01$Genename[sig01$padj<.05],sig02$Genename[sig02$padj<.05]))
  ## ranked by MDD
  sig <- left_join(sig01, sig02, by="Genename")
  sig <- sig[order(sig$padj.x),]
  degs[[area]] <- sig$Genename[which(sig$padj.x < 0.05 & sig$padj.y < 0.05)]
  print(area)
  print(sum(sig01$padj<.05))
  print(sum(sig02$padj<.05))
}

df <- data.frame(BA=c(9,11,24,25),MDD=c(367,3134,83,4065),PTSD=c(393,170,74,1), overlap=c(111,67,14,0))
df$percent.mdd <- (df$overlap/df$MDD * 100) %>% round(digits=3)
df$percent.ptsd <- (df$overlap/df$PTSD * 100) %>% round(digits=3)
write.csv(df, file="results/pc&degs/summary_mdd_ptsd.csv")






### WGCNA of miRNA ####
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


### for YI_UCI only
sel.ran <- sample(1:nrow(dat), 1000)
dat <- dat[sel.ran,]

## module construction
# R2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
# power <- min(which(R2>.8)) ## power=5
power <- 6
net <- blockwiseModules(datExpr = dat, power = power,
                        TOMType = "signed", minModuleSize = 10,
                        reassignThreshold = 0.2, mergeCutHeight = 0,
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
# datMeta1 <- data.frame(ptsd=as.numeric(datMeta$PrimaryDx=="PTSD"))
# datMeta1 <- data.frame(mdd=as.numeric(datMeta$PrimaryDx=="MDD"))
# datMeta1$mdd[datMeta$PrimaryDx=="PTSD"] <- NA
# datMeta1$ptsd[datMeta$PrimaryDx=="MDD"] <- NA
# datMeta1 <- cbind(datMeta1, datMeta0)
datMeta1 <- datMeta[c('Sex','AgeDeath','RIN','PMI')]
datMeta1$MDD <- as.numeric(datMeta$PrimaryDx=="MDD")
datMeta1$PTSD <- as.numeric(datMeta$PrimaryDx=="PTSD")
datMeta1$Sex <- as.numeric(datMeta1$Sex)
datMeta1$Race.AA <- as.numeric(datMeta$Race=="AA")
datMeta1$Race.EA <- as.numeric(datMeta$Race=="CAUC")
# moduleColors <- labels2colors(net$colors)
nGenes = ncol(dat);
MEs0 = moduleEigengenes(dat, net$colors)$eigengenes
MEs = orderMEs(MEs0)
nSamples <- dim(datMeta1)[1]
moduleTraitCor = cor(MEs, datMeta1, use = "p");# Pearson correlation
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# pdf('result_0703/wgcna_traitcor_dl_cp.pdf')
# pdf('result_0703/wgcna_traitcor_dl_cm.pdf')
# pdf('result_0703/wgcna_traitcor_sg_cp.pdf')
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
# dev.off()

nSamples = nrow(dat)
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
all <- data.frame(Geneid=datProbes2$Genename,
                  Protein=datProbes2$Geneid,
                  module=net$colors)
all <- cbind(all, trait=geneTraitSignificance)
mm <- sapply(1:dim(all)[1], 
             function(x){geneModuleMembership[x, match(all$module[x], modNames)]})
all$DxmoduleMembership <- mm

## save
save(all, datFPKM, datExpr2, datProbes2, datMeta, MEs, moduleTraitCor, moduleTraitPvalue,
     file = "result_2022/WGCNA_CP_ds_combatBatch.RData")

## follow up
lm.gs <- summary(lm(data=subset(all,module==1), trait.GS.PTSD ~ DxmoduleMembership))
plot(data=subset(all,module==1), trait.GS.PTSD ~ DxmoduleMembership, pch=16,
     main=paste0("R2=",lm.gs$adj.r.squared %>% round(digits=3)))
points(data=subset(all,Geneid %in% c("MIR6786","MIR589")), trait.GS.PTSD ~ DxmoduleMembership,pch=16,col="red")
# text(data=subset(all,Geneid %in% c("MIR6786","MIR589")), trait.GS.PTSD ~ DxmoduleMembership, Geneid)
abline(lm.gs$coefficients[,1], lty=2, col="red")
