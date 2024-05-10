rm(list=ls())
.rs.restartR()
set.seed(67)
gc()

library(Seurat)
library(dreamlet)
library(muscat)
library(tidyselect)
library(ExperimentHub)
library(zenith)
library(scater)
library(enrichR)
library(MAST)
library(NMF)
library(data.table)
library(EnhancedVolcano)
library(ggvenn)
library(BayesPrism)
library(DESeq2)
library(glue)
library(stringr)
library(stats)
library(dplyr)

#DE on the four major vascular cell types (i.e FB, SMC, Arterial EC, FBMC)
obj <- readRDS('Final_ax_object.RDS')
Idents(obj) <- obj$ax_status
markers.list <- c()
for(i in c('FB','SMC','FBMC')){
  print(i)
  subset <- subset(obj,subset = filtered_annotations == i)
  markers <- FindMarkers(
    subset, ident.1 = "Unruptured", ident.2 = "CTRL", 
    test.use = "MAST", verbose = T ,
    latent.vars = c("Sample",'Gender','Preparation',"Age"),
    max.cells.per.ident = 3000,
    random.seed = 67,
    densify = T)
  markers.list <- c(markers.list,list(markers)
  )
}
names(markers.list) <- c('FB','SMC','FBMC')

obj2 <- readRDS(file = "Final_EC_object.RDS")
Idents(obj2) <- obj2$ax_status
markers2 <- FindMarkers(
  obj2, ident.1 = "Unruptured", ident.2 = "CTRL", 
  test.use = "MAST", verbose = T ,
  latent.vars = c("Sample",'Gender','Preparation',"Age"), 
  max.cells.per.ident = 3000,
  random.seed = 67,
  densify = T
)

marker.list[['Art']] <- marker.list2
res.total <- c()
Total.cells <- c()
for(i in c('FB','SMC','Art','FBMC')) {
  print(i)
  celltype <- i
  scMAST_raw <- marker.list[[celltype]]
  res.filtered <- scMAST_raw[which(!(str_detect(rownames(scMAST_raw), "^ENSG")) &  !(str_detect(rownames(scMAST_raw), "^MT-")) & !(str_detect(rownames(scMAST_raw),"^RPS")) & !(str_detect(rownames(scMAST_raw), "^RPL")) ),]
  Up <- res.filtered %>% #up
    filter(p_val_adj < 0.01 & avg_log2FC > 1 & pct.1 > 0.1 & pct.2 > 0.1)
  Down <- res.filtered %>% #down
    filter(p_val_adj < 0.01 & avg_log2FC < -1 & pct.1 > 0.1 & pct.2 > 0.1)
  Total <- rbind(Up,Down) 
  Total$celltype <- i
  Total$gene <- rownames(Total)
  rownames(Total) <- NULL
  Total.cells <- rbind(Total.cells,Total)
  res.filtered$celltype <- i
  res.filtered$gene <- rownames(res.filtered)
  rownames(res.filtered) <- NULL
  res.total <- rbind(res.total,res.filtered)
}
res.total$sig <- 0.1
res.total$sig[which(res.total$p_val_adj < 0.01 & abs(res.total$avg_log2FC) > 1 & res.total$pct.1 > 0.1 & res.total$pct.2 > 0.1)] <- 1
res.total %<>% 
  mutate(sigxcelltype =  case_when(sig == 1 ~ celltype))
x_lim <- max(abs(-c(max(res.total$avg_log2FC),min(res.total$avg_log2FC))))

pruned.DEGs <- res.total %>% 
  filter(!is.na(sigxcelltype))

#Gene Burden
obj <- readRDS('Final_ax_object.RDS')
custom_palette <- unique(c('#0085ff',#FB
                           "#b9ceff",#SMC 
                           '#00c13c',#FBMC
                           "#225ba0"#EC
))
names(custom_palette) <- c('FB','SMC','FBMC','Art')
gene_burden <- c()
clusters <- c()
DEGs.clusters <- c()
DEGs.list <- c()
for(i in c('FB','SMC','FBMC')){
  Cells <- subset(obj,idents = i)
  if("CTRL" %in% unique(Cells@meta.data$ax_status) & "Unruptured" %in% unique(Cells@meta.data$ax_status)){
    if(table(Cells@meta.data$ax_status)["CTRL"]>250 & table(Cells@meta.data$ax_status)["Unruptured"]>250){
      for(j in 1:10){
        print(paste0(i,"_",j))
        Idents(Cells) <- Cells@meta.data$ax_status
        Cells.subset1 <- subset(Cells,idents = c("CTRL"))
        Cells.subset2 <- subset(Cells,idents = c("Unruptured"))
        set.seed(j)
        index1 <- sample(1:ncol(Cells.subset1),3000,replace=F)
        index2 <- sample(1:ncol(Cells.subset2),200,replace=F)
        cellid1 <- colnames(Cells.subset1)[index1]
        cellid2 <- colnames(Cells.subset2)[index2]
        Cells.subset <- subset(Cells,cells = c(cellid1,cellid2))
        Cells.DEGs <- FindMarkers(
          Cells.subset, ident.1 = "Unruptured", ident.2 = "CTRL", 
          test.use = "MAST", verbose = T ,
          latent.vars = c("Sample",'Gender','Preparation',"Age")
        )
        Cells.DEGs %<>%
          filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1 & pct.1 > 0.1 & pct.2 > 0.1)
        Cells.DEGs <- Cells.DEGs[which(!(str_detect(rownames(Cells.DEGs), "^ENSG")) &  !(str_detect(rownames(Cells.DEGs), "^MT-")) & !(str_detect(rownames(Cells.DEGs),"^RPS")) & !(str_detect(rownames(Cells.DEGs), "^RPL")) ),]
        gene_burden <- c(gene_burden,nrow(Cells.DEGs)) 
        clusters <- c(clusters,i)
        DEGs.clusters <- c(DEGs.clusters,paste0(i,"_",j))
        DEGs.list <- c(DEGs.list,list(Cells.DEGs))
        print(nrow(Cells.DEGs))
      }
    }else{print(paste0("Skipping_",i))}
  }else{print(paste0("Skipping_",i))}
}
gene_burden.merged <- as_tibble(cbind(gene_burden,clusters))
gene_burden.merged$clusters <- as.factor(gene_burden.merged$clusters)
ordered <- gene_burden.merged %>% 
  group_by(clusters) %>% 
  summarise(mean = mean(as.numeric(gene_burden))) %>% 
  ungroup() %>% 
  arrange(by_group = mean)
gene_burden.merged$clusters <- fct_relevel(gene_burden.merged$clusters,rev(as.character(ordered$clusters)))
names(DEGs.list) <- DEGs.clusters

DEGs.ledger <- c()
for(i in 1:10){
  DEGs.ledger <- c(DEGs.ledger,rownames(as.data.frame(DEGs.list[1])))
}
DEGs.ledger <- unique(DEGs.ledger)

# Arterialized
obj <- readRDS('Final_EC_object.RDS')
table(obj$ax_status)
gene_burden <- c()
clusters <- c()
DEGs.clusters <- c()
DEGs.list <- c()
for(i in c('Art')){
  Cells <- obj
  if("CTRL" %in% unique(Cells@meta.data$ax_status) & "Unruptured" %in% unique(Cells@meta.data$ax_status)){
    if(table(Cells@meta.data$ax_status)["CTRL"]>250 & table(Cells@meta.data$ax_status)["Unruptured"]>250){
      for(j in 1:10){
        print(paste0(i,"_",j))
        Idents(Cells) <- Cells@meta.data$ax_status
        Cells.subset1 <- subset(Cells,idents = c("CTRL"))
        Cells.subset2 <- subset(Cells,idents = c("Unruptured"))
        set.seed(j)
        index1 <- sample(1:ncol(Cells.subset1),3000,replace=F)
        index2 <- sample(1:ncol(Cells.subset2),200,replace=F)
        cellid1 <- colnames(Cells.subset1)[index1]
        cellid2 <- colnames(Cells.subset2)[index2]
        Cells.subset <- subset(Cells,cells = c(cellid1,cellid2))
        Cells.DEGs <- FindMarkers(
          Cells.subset, ident.1 = "Unruptured", ident.2 = "CTRL", 
          test.use = "MAST", verbose = T ,
          latent.vars = c("Sample",'Gender','Preparation',"Age")
        )
        Cells.DEGs %<>%
          filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1 & pct.1 > 0.1 & pct.2 > 0.1)
        Cells.DEGs <- Cells.DEGs[which(!(str_detect(rownames(Cells.DEGs), "^ENSG")) &  !(str_detect(rownames(Cells.DEGs), "^MT-")) & !(str_detect(rownames(Cells.DEGs),"^RPS")) & !(str_detect(rownames(Cells.DEGs), "^RPL")) ),]
        gene_burden <- c(gene_burden,nrow(Cells.DEGs)) 
        clusters <- c(clusters,i)
        DEGs.clusters <- c(DEGs.clusters,paste0(i,"_",j))
        DEGs.list <- c(DEGs.list,list(Cells.DEGs))
        print(nrow(Cells.DEGs))
      }
    }else{print(paste0("Skipping_",i))}
  }else{print(paste0("Skipping_",i))}
}
gene_burden.merged2 <- as_tibble(cbind(gene_burden,clusters))
gene_burden.merged2$clusters <- as.factor(gene_burden.merged2$clusters)
ordered <- gene_burden.merged2 %>% 
  group_by(clusters) %>% 
  summarise(mean = mean(as.numeric(gene_burden))) %>% 
  ungroup() %>% 
  arrange(by_group = mean)
gene_burden.merged2$clusters <- fct_relevel(gene_burden.merged2$clusters,rev(as.character(ordered$clusters)))
names(DEGs.list) <- DEGs.clusters

# Getting Gene Burdens
all.gene_burdens <- rbind(gene_burden.merged1,gene_burden.merged2)
all.gene_burdens %<>% 
  group_by(clusters) %>% 
  summarise(mean = mean(as.numeric(gene_burden)),
            gene_burden = gene_burden
  ) %>% 
  ungroup()
aov.results <- aov(gene_burden ~ clusters,data = all.gene_burdens) #stats
summary(aov.results)
tukey.two.way <- TukeyHSD(aov.results) #stats
tukey.two.way


# DE on rupt+ aMacs
obj <- readRDS('Final_ax_object.RDS')
obj$Scissor3 <- obj$Scissor
Idents(obj) <- obj$Scissor3
markers.list <- c()
for(i in c('aMac','PVM')){
  print(i)
  subset <- subset(obj,subset = filtered_annotations == i)
  markers <- FindMarkers(
    subset, ident.1 = "Rupture", ident.2 = "BG", 
    test.use = "MAST", verbose = T ,
    latent.vars = c("Sample",'Gender','Preparation',"Age"),
    max.cells.per.ident = 3000,
    random.seed = 67,
    densify = T)
  markers.list <- c(markers.list,list(markers)) 
}
names(markers.list) <- c('aMac','PVM')

aMac <- markers.list$aMac %>% 
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1 & pct.1 > 0.1 & pct.2 > 0.1)
PVM <- markers.list$PVM %>% 
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1 & pct.1 > 0.1 & pct.2 > 0.1)