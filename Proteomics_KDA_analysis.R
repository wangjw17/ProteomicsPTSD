##### This script is written for calculating KDA for proteomics data to avoid cache problem
##### Jiawei Wang, Dec 27, 2021
### KDA set up ####
library('grDevices')
library('WGCNA')
library('bnlearn')
library('KDA')
library('qgraph')
library('plyr')
library('pracma')
library('ggplot2')
library('YaleToolkit')
library('igraph')
library('circlize')
library('tidyr')
library('dplyr')
library('grDevices')
library('ggplotify')

### directories and locations
rm(list=ls())
# args <- c("~/Documents/Projects/PTSD/proteomics/DIA_new.RData", ## DATA FILE
#           "~/Documents/Projects/PTSD/proteomics/",  ## DEG_FILE_PATH
#           "result_0518/DE_MDD_PTSD_dl.csv",  ## DEG_FILE_PATTERN
#           "~/Documents/Projects/PTSD/proteomics/",  ## WGCNA_FILE_PATH
#           "result_0703/WGCNA_CP_dl.RData",  ## WGCNA_FILE_PATTERN
#           "result_0128/",  ## OUTPUT_FILE_PATH
#           "0208") ## OUTFILE_SUFFIX
# Figures_path <- args[6]
# outfile_suffix <- args[7]

### parameteres needed
## degs to all0102
# load('DIA_new.RData')
load('DIA_new1.RData')
fpkmdat <- as.data.frame(datExpr)
demos <- datMeta
fpkmdat$Genename <- datProbes$GN
# dep <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
dep <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
thres_padj <- 0.05 ## args[8]
all0102 <- data.frame(Geneid = dep$GN, 
                      MDD.A9 = as.numeric(dep$MDD.P.Value < thres_padj) * sign(dep$MDD.logFC),
                      PTSD.A9 = as.numeric(dep$PTSD.P.Value < thres_padj) * sign(dep$PTSD.logFC))
names(all0102)[1] <- "Genename"
group <- "A9"


## enrichment of DEPs
# sig.p <- read.csv('result_0518/DE_MDD_PTSD_dl.csv')
sig.p <- read.csv('result_0518/DE_MDD_PTSD_dl_incA_new1.csv')
dep.m <- sig.p$GN[sig.p$MDD.P.Value<.05]
dep.p <- sig.p$GN[sig.p$PTSD.P.Value<.05]

# sig.p <- read.csv('result_0128/DE_PvCM_dlsg.csv')
# dep.m <- sig.p$GN[sig.p$dlPFC.P.Value<.05]
# dep.p <- sig.p$GN[sig.p$sgPFC.P.Value<.05]

## load WGCNA results
# load('result_0703/WGCNA_CM_dl.RData'); names(all)[1] <- "Genename" ##Change "Geneid" to "Genename";
load('result_0703/WGCNA_CP_dl.RData')

colnames(moduleTraitPvalue)[1] <- "ptsd"

N <- (all$Genename %>% unique) %in% dep.m %>% sum
modules <- levels(factor(all$module))
MM <- nrow(all)
allm <- data.frame(modules=modules, enrich=NA, pval=NA)
for (m in modules){
  M <- sum(all$module==m)
  n <- (dep.m %in% all$Genename[all$module==m]) %>% sum
  allm$enrich[allm$modules==m] <- n*MM/(M*N)
  allm$pval[allm$modules==m] <- fisher.test(table(data.frame(
    A=all$module==m, B=all$Genename %in% dep.m)))$p.value
}
allm1 <- allm; allm1$dx = "MDD"
N <- (all$Genename %>% unique) %in% dep.p %>% sum
modules <- levels(factor(all$module))
MM <- nrow(all)
allm <- data.frame(modules=modules, enrich=NA, pval=NA)
for (m in modules){
  M <- sum(all$module==m)
  n <- (dep.p %in% all$Genename[all$module==m]) %>% sum
  allm$enrich[allm$modules==m] <- n*MM/(M*N)
  allm$pval[allm$modules==m] <- fisher.test(table(data.frame(
    A=all$module==m, B=all$Genename %in% dep.p)))$p.value
}
allm2 <- allm; allm2$dx = "PTSD"
allm <- rbind(allm1, allm2)
modules.enrichment <- data.frame(MDD.A9 = allm$pval[allm$dx=="MDD"],PTSD.A9 = allm$pval[allm$dx=="PTSD"])
rownames(modules.enrichment) <- allm$modules[allm$dx=="PTSD"]

## WGCNA results
# dend <- denro.row
minfo <- as.data.frame(all[, c('Geneid', 'module', 'Genename')], stringsAsFactors = F)
colnames(minfo) <- c('Geneid', 'module', 'Genename')
minfo <- join(minfo, all0102, by = 'Genename', 'left', 'first')
# Print the number of DEGs 
print(sprintf('For this dataset, there are %d genes differentially expressed for MDD and %d genes for PTSD', 
              sum(abs(minfo[, paste0('MDD.', group)]), na.rm = T), sum(abs(minfo[, paste0('PTSD.', group)]), na.rm = T)))

rownames(moduleTraitPvalue) <- gsub('ME', '', rownames(moduleTraitPvalue))
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
moduleTraitPvalue$Dx.PTSD.fdr <- p.adjust(moduleTraitPvalue$ptsd, 'fdr')
modules_csea_fdr <- to_plot
modules <- names(table(minfo$module))
moduleTraitPvalue <- moduleTraitPvalue[modules, ]
modules_csea_fdr <- modules_csea_fdr[modules, ]


### KDA analysis simplified ####
colors <- matrix(c('brown1', 'white', 'brown4',               # Color for 4
                   'cyan1', 'white', 'darkblue',              # Color for 9
                   'darkolivegreen1', 'white', 'forestgreen', # Color for 11
                   'gold', 'white', 'gold4',                  # Color for 24
                   'azure3', 'white', 'azure4' ),             # Color for 25
                 ncol= 3, byrow = T) 


# pdf('result_0128/kda/KDA_CM_dl_new1.pdf')
pdf('result_0128/kda/KDA_CP_dl_new1.pdf')
for(m in modules){ ## CM/CP_all
  DEm <- minfo[which(minfo$module == m & minfo[, paste0('PTSD.', group)]), 'Genename', drop=F] ## gene names
  datExpr <- join(minfo[minfo$module == m & !is.na(minfo[, paste0('PTSD.', group)]), c('Geneid', 'Genename'), drop = F], 
                  fpkmdat[, c('Genename', demos$Sample[which(demos$Region == "dlPFC")])], 
                  by = 'Genename', 'left', 'first')
  head(datExpr)
  datET <- collapseRows(datET = datExpr[, -grep('Gene', colnames(datExpr))], rowGroup = datExpr[, 'Genename'], 
                        rowID = rownames(datExpr), connectivityBasedCollapsing = T)$datETcollapsed
  ## reduce genename duplicity; #gene x #sample
  plot_genes <- rownames(datET) 
  print(paste0('Module ', m, ' has ', dim(datET)[1], ' genes.'))
  
  # Infer the graph using aracne
  aracne.results <- aracne(as.data.frame(t(datET)))
  connection_m <- aracne.results$arcs
  connection_m <- connection_m[connection_m[, 2] < connection_m[, 1], ,drop=F]
  
  # Run key driver analysis based on the graph
  key.drivers <- list()
  if(dim(datExpr)[1] >= 10){
    KDA_results <- KDA::keyDriverAnalysis(connection_m, signature = as.character(DEm[,1]), directed = F, nlayer_search=3, expanded_network_as_signature=T)
    if(is.null(KDA_results)){
      cols <- c('keydrivers','is_signature','hits','downstream','signature_in_subnetwork','subnetwork_size','signature_in_network','network_size','signature','optimal_layer','fold_change_whole','pvalue_whole','fold_change_subnet','pvalue_subnet','pvalue_corrected_subnet','keydriver','Module')
      key.driver <- as.data.frame(matrix(NA, nrow = 0, ncol = length(cols)))
      colnames(key.driver) <- cols
    }
    else{
      key.driver <- as.data.frame(KDA_results$keydrivers, stringsAsFactors = F)
      key.driver$Module <- m
      key.driver <- key.driver[key.driver$is_signature == 1, ]
      key.drivers[[length(key.drivers)+1]] <- key.driver
    }
  }
  if(dim(datET)[1] > 50){
    DEonly <- apply(aracne.results$arcs, 1, function(x) (x[1] %in% DEm[, 1]) | (x[2] %in% DEm[, 1]))
    connection_m <- aracne.results$arcs[DEonly,,drop=F]
    connection_m <- connection_m[connection_m[, 2] < connection_m[, 1], ,drop=F]
    plot_genes <- unique(c(aracne.results$arcs[DEonly, 'from'], aracne.results$arcs[DEonly, 'to']))
  }
  if(length(plot_genes) > 100){
    DEonly <- apply(connection_m, 1, function(x) (x[1] %in% c(key.driver$keydrivers, DEm[, 1])) & (x[2] %in% c(key.driver$keydrivers, DEm[, 1])))
    connection_m <- connection_m[DEonly, ,drop=F]
    plot_genes <- unique(as.character(connection_m))
  }
  
  # Get the annotation
  csea_anno <- paste0('CSEA: ', paste0(colnames(modules_csea_fdr)[modules_csea_fdr[m, ] < 0.05], collapse = ', '))
  trait_anno <- paste0('Trait: ', c('PTSD')[moduleTraitPvalue[m,'ptsd'] < 0.05]) ## for PTSD
  # trait_anno <- paste0('Trait: ', c('MDD')[moduleTraitPvalue[m,'ptsd'] < 0.05]) ## for MDD
  # Plot the network
  graph_m <- graph_from_edgelist(connection_m, directed = F)
  vertex.frame.color <- 'black'
  v.sizes <- linspace(x1 = 4, x2 = 20, 40)
  
  vertex.frame.color <- lapply(V(graph_m)$name, function(x) {
    if(x %in% key.driver$keydrivers){
      return('magenta')
    }
    else return('black')
  })
  vertex.frame.color <- do.call(c, vertex.frame.color)
  V(graph_m)$label.cex = .6
  
  ## for single region
  pie.values <- lapply(V(graph_m)$name, function(x) rep(1, length(names(which(table(demos$Region[demos$Region=="dlPFC"]) > 0)))))
  pie.color <- lapply(V(graph_m)$name, function(x) {
    tmp <- as.numeric(minfo[minfo$Genename == x, grep('PTSD.[A-Z]+[0-9]+', colnames(minfo))] + 2)
    return(colors[as.matrix(cbind(1:length(tmp), tmp))])
  })
  # ## for multiple regions/groups
  # pie.values <- lapply(V(graph_m)$name, function(x) rep(1, length(names(which(table(demos$Region) > 0)))))
  # pie.color <- lapply(V(graph_m)$name, function(x) {
  #   tmp <- as.numeric(as.matrix(minfo[minfo$Genename == x, grep('A9', colnames(minfo))] + 2)) ## some genes with multiple entries
  #   return(colors[as.matrix(cbind(1:length(tmp), tmp))])
  # })
  
  pvalue <- sprintf("%.2e", modules.enrichment[m, 'PTSD.A9']) ## PTSD
  # pvalue <- sprintf("%.2e", modules.enrichment[m, 'MDD.A9']) ## MDD
  plot(graph_m, 
       vertex.label.dist=igraph::degree(graph_m) * 0.05 + 0.5,
       vertex.label.degree=pi/2,
       vertex.shape='pie', 
       vertex.frame.color = vertex.frame.color,
       vertex.pie = pie.values,
       vertex.pie.color=pie.color,
       vertex.size = v.sizes[igraph::degree(graph_m)]) +
    title(main = paste0('Module ', m, " of PTSD (p-value=", pvalue, ") \n", csea_anno, "\n", trait_anno)) ## for PTSD
    # title(main = paste0('Module ', m, " of MDD (p-value=", pvalue, ") \n", csea_anno, "\n", trait_anno)) ## for MDD
}
dev.off()

