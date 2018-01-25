library(DESeq2)
library(pamr)
library(biomaRt)
mart=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",mart=mart)

load('~/Honours/Stage-Prediction-of-Cancer/papillary/environment/dds.RData')
load('~/Honours/Stage-Prediction-of-Cancer/ccrcc/environment/kirc_data.RData')
met.genes.df <- read.csv('metabolic_genes.csv')
met.genes <- as.character(met.genes.df$GENE.ID.1)
genes.entrez = getBM(attributes = c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = g, mart = ensembl)


rownames(dds) <- remove.dots(rownames(dds))
length(intersect(rownames(prcc.data), met.genes)) == length(met.genes)
length(intersect(rownames(data), met.genes)) == length(met.genes)

met.df.prcc <- assay(dds)[met.genes,]


remove.dots <- function(ens.ids.all)
{
  ###ens.ids.all <- gets the ids returned from get.genes.files
  
  ##The ens ids contain symbols after dots making them as invalid ensembl ids for using for enrichment
  ##analysis, so stripping the same ids removing the unwanted things after dot
  
  g = sapply(ens.ids.all, function(x) 
  {
    unlist(strsplit(x, split = '.', fixed = T))[1]
  }) ##removing the symbols after .
  return(g)
}

prcc.matched.data <- prcc.data[,match(colnames(dds), colnames(prcc.data))]
prcc.matched.data <- prcc.matched.data[-which(rowSums(assay(prcc.matched.data)) < 10),]

prcc.matched.data.met <- prcc.matched.data[intersect(met.genes, rownames(prcc.matched.data)),]
dds.obj <- DESeqDataSetFromMatrix(assay(prcc.matched.data.met),
                colData = colData(prcc.matched.data.met), design = ~shortLetterCode)
dds.obj.ent <- DESeqDataSetFromMatrix(assay(prcc.matched.data),
                colData = colData(prcc.matched.data), design = ~shortLetterCode)
dds.obj <- DESeq(dds.obj, parallel = T) 
dds.obj.ent <- DESeq(dds.obj.ent, parallel = T)

res <- results(dds.obj, contrast = c('shortLetterCode', 'TP', 'NT'), parallel = T)
res.ent <- results(dds.obj.ent, contrast = c('shortLetterCode', 'TP', 'NT'), parallel = T)
summary(res)
summary(res.ent)
g <- get.genes(res.prcc[[1]], 2, 0.05, 0.05)
g.ent <- intersect(get.genes(res.prcc[[2]], 2, 0.05, 0.05), met.genes)

prcc.pat <- get.matched.ind(colData(prcc.data))
ccrcc.pat <- get.matched.ind(colData(data))

res.prcc <- get.deseq2.proc(prcc.data, prcc.pat, met.genes)
res.ccrcc <- get.deseq2.proc(data, ccrcc.pat, met.genes)

res.sam.prcc <- get.sam.res(prcc.data, prcc.pat, met.genes)
res.sam.ccrcc <- get.sam.res(data, ccrcc.pat, met.genes)

g.sam.prcc <- sapply(get.sam.genes(res.sam.prcc, list(2,3,4,5)), function(g) intersect(g, met.genes))
g.sam.ccrcc <- sapply(get.sam.genes(res.sam.ccrcc, list(2,3,4,5)), function(g) intersect(g, met.genes))
g.sam.chcc <- sapply(get.sam.genes(res.sam.ch, list(2,3,4,5)), function(g) intersect(g, met.genes))

g.deseq.prcc <- lapply(c(2,3,4,5), function(x){intersect(get.deseq2.genes(res.prcc[[2]], x, 0.05, 0.05), met.genes)})
names(g.deseq.prcc) <- c('2 fold', '3 fold', '4 fold', '5 fold')
g.deseq.ccrcc <- lapply(c(2,3,4,5), function(x){intersect(get.deseq2.genes(res.ccrcc[[2]], x, 0.05, 0.05), met.genes)})
names(g.deseq.ccrcc) <- c('2 fold', '3 fold', '4 fold', '5 fold')
g.deseq.chr <- lapply(c(2,3,4,5), function(x){intersect(get.deseq2.genes(res.deseq.ch[[2]], x, 0.05, 0.05), met.genes)})
names(g.deseq.chr) <- c('2 fold', '3 fold', '4 fold', '5 fold')

g.deseq.prcc.met <- lapply(c(2,3,4,5), function(x){get.deseq2.genes(res.prcc[[1]], x, 0.05, 0.05)})
names(g.deseq.prcc.met) <- c('2 fold', '3 fold', '4 fold', '5 fold')

library(pheatmap)
ann.col.df <- data.frame(type=colData(prcc.data)$shortLetterCode[unlist(prcc.pat)], 
                         row.names = colnames(prcc.data)[unlist(prcc.pat)])
pheatmap(assay(prcc.data)[g[1:10], unlist(prcc.pat)], cluster_rows = F, 
         cluster_cols = T, annotation_col = ann.col.df, show_colnames = F )
