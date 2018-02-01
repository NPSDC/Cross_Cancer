source('download.R')

met.genes.df <- read.csv('metabolic_genes.csv')
met.genes <- as.character(met.genes.df$GENE.ID.1)

data.canc <- download.data(project.id = 'TCGA-KIRP', directory = '~/GDC_Data/')
match.ind <- get.matched.ind(col.data = colData(data.canc))
res.deseq <- get.deseq2.proc(sum.exp = data.canc, match.ind = match.ind, genes.proc = met.genes)
res.sam <- get.sam.res(sum.exp = data.canc, match.ind = match.ind, genes.proc = met.genes)

g.sam <- get.sam.genes(res.sam.df = res.sam, folds = list(2,3,4,5))
g.deseq <- get.deseq2.genes(res = res.deseq, folds = list(2,3,4,5))
