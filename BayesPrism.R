rm(list=ls())
.rs.restartR()
set.seed(67)
gc()

library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(dplyr)
library(tidyr)
library(irlba)
library(stringr)
library(EnhancedVolcano)
library(gprofiler2)
library(forcats)
library(speckle)
library(limma)
library(UCell)
library(BayesPrism)
library(DESeq2)
library(variancePartition)
library(VennDiagram)
library(monocle3)
library(SeuratWrappers)
library(richR)
library(MAST)
library(dittoSeq)
library(viridis)
library(infer)
library(RColorBrewer)
library(ArchR)
library(enrichR)
library(EnhancedVolcano)
library(Scissor)
library(ggdendro)
library(pheatmap)
library(magrittr)
library(ggvenn)

# Deconvolutions on 
bk.dat <- read.csv("Final_bulk_counts_matrix.csv",row.names = NULL)
duplicates <- bk.dat$row.names[duplicated(bk.dat$row.names)]
bk.dat <- bk.dat[-which(bk.dat$row.names %in% duplicates),]
rownames(bk.dat) <- bk.dat$row.names
bk.dat <- t(bk.dat[,-1])

obj <- readRDS(file='Final_ax_object.RDS')
sc.dat <- obj
sc.dat <- t(as.data.frame(sc.dat@assays[["SCT"]]@counts))
cell.state.labels <- as.character(obj$filtered_annotations)
cell.type.labels <- cell.state.labels
cell.type.labels[which(cell.type.labels %in% c('PVM','Mono','aMac','cDC','MG'))] <- 'Myeloid'

table(cbind.data.frame(cell.state.labels, cell.type.labels))

#Starting Bayes Prism
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell type correlation",
              cexRow=0.7, cexCol=0.7,
              margins=c(4,4))

sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=cell.state.labels,
  species="hs", 
  return.raw=TRUE
)

head(sc.stat)  

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=sc.dat,
  cell.type.labels=cell.state.labels,
  species="hs",
  return.raw=TRUE
)

head(bk.stat)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
)

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=120)
saveRDS(bp.res,file = "Final_BayesPrism_Model.RDS")





