library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(samr)
source('sam_func_cop.R')
source('samr.morefuns.R')

download.data <- function(project.id, directory)
{
  down.dir <- paste(directory, project.id)
  if(project.id %in% getGDCprojects()$project_id)
  {
    q <- GDCquery(project = project.id, data.category = 'Transcriptome Profiling', 
                  data.type = 'Gene Expression Quantification', 
                  workflow.type = 'HTSeq - Counts')
    if(dir.exists(directory))
      dir.create(paste(directory, project.id))
    else
    {
      print("Invalid directory path")
      return(0)
    }
    GDCdownload(query = q, directory = '~/PRCC/')
    data <- GDCprepare(q, directory = '~/PRCC')
    return(data) 
  }
  else
    print('Invalid Project')
}

get.matched.ind <- function(col.data)
{
  matched.nor.ind <- c()
  matched.tum.ind <- c()
  dup.ind <- which(duplicated(col.data$patient))
  tum.ind <- which(col.data$shortLetterCode == 'TP')
  nor.ind <- which(col.data$shortLetterCode == 'NT')
  for(i in seq_along(dup.ind))
  {
    ind.match <- which(col.data$patient == col.data$patient[dup.ind[i]])
    tum.int <- intersect(tum.ind, ind.match)
    nor.int <- intersect(nor.ind, ind.match)
    if(length(nor.int) == 0 | length(tum.int) == 0)
      cat("Matched normal tumor ID doesn't exist", ind.match)
    else
    {
      matched.nor.ind <- c(matched.nor.ind, nor.int)
      matched.tum.ind <- c(matched.tum.ind, tum.int)
    }
  }
  return(list(norm.ind = matched.nor.ind, tumor.ind = matched.tum.ind))
}

get.deseq2.proc <- function(sum.exp, match.ind, genes.proc)
{
  library("BiocParallel")
  register(MulticoreParam(4))
  data.mat <- assay(sum.exp)
  ex.count <- which(rowSums(data.mat) < 10)
  if(length(ex.count) > 0)
    data.mat <- data.mat[-ex.count,]
  type = rep(c('T'), ncol(data.mat))
  type[c(match.ind$tumor.ind)] = 'MT'
  type[c(match.ind$norm.ind)] = 'N'
  
  col.data <- colData(sum.exp)
  col.data$type = type
  
  #data.mat.genes.proc <- data.mat[intersect(genes.proc, rownames(data.mat)),]
  #ex.count <- which(rowSums(data.mat.genes.proc) < 10)
  #if(length(ex.count) > 0)
   # data.mat.genes.proc <- data.mat.genes.proc[-ex.count,]
  #data.mat <- data.mat[intersect(met.genes, rownames(data.mat)),]
  #dds.obj.genes.proc <- DESeqDataSetFromMatrix(data.mat.genes.proc, colData = col.data,
   #                                  design = ~type)
  dds.obj.ent <- DESeqDataSetFromMatrix(data.mat, colData = col.data, design = ~type)
                                    
  #dds.obj.genes.proc <- DESeq(dds.obj.genes.proc, parallel = T) 
  dds.obj.ent <- DESeq(dds.obj.ent, parallel = T)
  
  #res.genes.proc <- results(dds.obj.genes.proc, contrast = c('type', 'MT', 'N'), 
  #                          parallel = T)
  res.ent <- results(dds.obj.ent, contrast = c('type', 'MT', 'N'), parallel = T)
  res.ent <- res.ent[intersect(rownames(res.ent), genes.proc),]
  #return(list(res.genes.proc, res.ent))
  return(res.ent)
}

get.sam.res <- function(sum.exp, match.ind, genes.proc)
{
  type = rep(c(1), length(unlist(match.ind)))
  type[1:length(match.ind$norm.ind)] = 1
  type[(1+length(match.ind$norm.ind)):length(unlist(match.ind))] = 2
  
  sum.exp.matched <- sum.exp[,unlist(match.ind)]
  data.mat <- assay(sum.exp.matched)
  ex.count <- which(rowSums(data.mat) < 10)
  if(length(ex.count) > 0)
    data.mat <- data.mat[-ex.count,]

  res.sam <- SAMseq(x = data.mat, y = type, resp.type = 'Two class unpaired', 
                    genenames = rownames(data.mat))
  res.sam.df <- get.sam.df(res.sam)
  return(res.sam.df[intersect(rownames(res.sam.df), genes.proc),])
}

get.sam.df <- function(sam.obj)
{
  
    res.df <- c()
    if(sam.obj$siggenes.table$ngenes.up == 0)
      res.df <- sam.obj$siggenes.table$genes.lo
    else if(sam.obj$siggenes.table$ngenes.lo == 0)
      res.df <- sam.obj$siggenes.table$genes.up
    else
      res.df <- rbind(sam.obj$siggenes.table$genes.up, sam.obj$siggenes.table$genes.lo)
    res.df <- data.frame(GeneID = as.character(res.df[,1]), Gene.Name = as.character(res.df[,2]), 
               Score = as.numeric(res.df[,3]), Fold.Change = as.numeric(res.df[,4]),
               q.val = as.numeric(res.df[,5]), stringsAsFactors = F)
    rownames(res.df) = res.df$GeneID
  
  return(res.df)
}

get.sam.genes <- function(res.sam.df, folds)
{
  req.genes <- lapply(folds, function(fold)
  {
      genes <- res.sam.df$GeneID[which(abs(log2(res.sam.df$Fold.Change)) > fold & 
                                                   res.sam.df$q.val < 0.05)]
    
    
  })
  names(req.genes) <- paste0(folds, ' fold')
  return(req.genes)
}

get.deseq2.genes <- function(res, folds)
{
  req.genes <- lapply(folds, function(fold)  return(rownames(res)[abs(res[,2]) > fold & res[,6] < 0.05 & res[,5] < 0.05]))
  names(req.genes) <- paste0(folds, ' fold')
  return(req.genes)
}

get.gene.symbol <- function(genes, met.df)
{
  g.ind <- match(genes, met.df[,1])
  return(as.character(met.df[,4][g.ind]))
}