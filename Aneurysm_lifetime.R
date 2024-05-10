rm(list=ls())
.rs.restartR()
set.seed(67)
gc()

library(Seurat)
library(dplyr)
library(monocle3)
library(viridis)
library(ggplot2)
library(BayesPrism)
library(DESeq2)
library(magrittr)
library(stringr)
library(enrichR)
library(tidyr)
library(harmony)
library(UCell)
library(ggpubr)

bk.dat <- read.csv("Final_bulk_counts_matrix.csv",row.names = NULL)
bk.meta <- read.csv("Final_geo_aneurysm_meta_data.csv")
duplicates <- bk.dat$row.names[duplicated(bk.dat$row.names)]
bk.dat <- bk.dat[-which(bk.dat$row.names %in% duplicates),]
rownames(bk.dat) <- bk.dat$row.names
bk.dat <- bk.dat[,-1]
rownames(bk.meta) <- bk.meta$Run
bk.meta$rupture.status[which(is.na(bk.meta$rupture.status))] <- 'CTRL'

bk.seurat <- CreateSeuratObject(bk.dat,meta.data = bk.meta)
bk.seurat <- SCTransform(bk.seurat, verbose = T)
bk.seurat <- RunPCA(bk.seurat, verbose = T)
bk.seurat <- RunUMAP(bk.seurat, dims = 1:30, verbose = T,min.dist = 0.01) #Used for the bulk only
DimPlot(bk.seurat,group.by = 'rupture.status',pt.size = 5)

cds <- SeuratWrappers::as.cell_data_set(bk.seurat)
cds <- cluster_cells(cds,partition_qval = 0.05,random_seed = 67)
cds <- learn_graph(cds,learn_graph_control = list(minimal_branch_len = 8))
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE,cell_size = 5,color_cells_by = 'rupture.status',group_label_size = 5)
max <- which.max(unlist(FetchData(bk.seurat, "AVP")))
max <- colnames(bk.seurat)[max]
cds <- order_cells(cds, root_cells = max)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE,show_trajectory_graph = F,rasterize = F,cell_size = 5) + 
  scale_color_viridis()

# Pull the dx pseudotime
traj.plot <- plot_cells(cds, color_cells_by = "pseudotime")
point.data <- ggplot_build(traj.plot)[["plot"]][["data"]] %>% 
  as_tibble() %>% 
  select(Run,cell_color) %>% 
  arrange(cell_color)

#Creating the pseudotimed seurat object
bk.seurat@meta.data %<>% 
  left_join(point.data,by = 'Run')
rownames(bk.seurat@meta.data) <- bk.seurat$Run
bk.seurat$cell_color2 <- -as.numeric(bk.seurat$cell_color)
bk.seurat$cell_color2 <- bk.seurat$cell_color2 + max(bk.seurat$cell_color)
bk.seurat <- subset(bk.seurat,subset = tissue.status != 'Unknown intracranial aneurysm')
DimPlot(bk.seurat,group.by = 'tissue.status',pt.size = 5) + NoLegend() + xlim(-5,5) + ylim(-3,3)
FeaturePlot(bk.seurat,features = 'cell_color2',pt.size = 5) + NoLegend() + xlim(-5,5) + ylim(-3,3) +
  scale_color_viridis()  
saveRDS(bk.seurat,'Final_bulk_dx_pseudotime_obj.RDS')



# Looking at celltype changes across dx pseudotime
DEGs <- read.csv('Final_annotated_markers.csv',row.names = 1)
DEGs %<>% 
  as_tibble() %>% 
  group_by(cluster) %>% 
  slice_max(n = 25,order_by = weighted_specificity) %>% 
  ungroup() %>%
  select(cluster,gene) %>%
  pivot_wider(names_from = cluster,values_from = gene,values_fn = list)
DEGs <- as.list(DEGs)

bk.seurat <-  readRDS('Final_bulk_dx_pseudotime_obj.RDS')
bk.seurat <- AddModuleScore_UCell(bk.seurat,features=DEGs)
bk.dat <- as.matrix(bk.seurat@assays$SCT$data)
read_depth <- as.data.frame(colSums(bk.dat))
bk.meta <- bk.seurat@meta.data
bk.meta %<>% 
  as_tibble() %>% 
  mutate(encoded = case_when(rupture.status == 'Ruptured' ~ 2,
                             rupture.status == 'Unruptured' ~ 1,
                             condition == 'control' ~ 0))
gene_signature <- as.data.frame(cbind(bk.meta$Run,as.numeric(bk.meta$aFB1_UCell),as.numeric(bk.meta$aMac_UCell),as.numeric(bk.meta$SMC_UCell)))
colnames(gene_signature) <- c('Run','aFB1','aMac','SMC')
gene_signature_values <- gene_signature[,2:ncol(gene_signature)]
for(i in 1:ncol(gene_signature_values)){
  seed_idx <- which(gene_signature$Run == seed_run)
  setpoint <- gene_signature_values[seed_idx,i]
  gene_signature_values[,i] <- as.numeric(gene_signature_values[,i]) - as.numeric(setpoint)
}
gene_signature <- cbind(gene_signature$Run,gene_signature_values)
colnames(gene_signature)[1] <- 'Run'
gene_signature %<>% 
  as_tibble() %>% 
  pivot_longer(cols = c('aFB1','aMac','SMC'),names_to = 'celltype')
bk.meta %>% 
  left_join(gene_signature,by = 'Run') %>% 
  ggplot(aes(cell_color2,as.numeric(value))) +
  geom_point(aes(color = celltype)) + 
  geom_smooth(aes(color = celltype),method='loess',se = F) + 
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_color_manual(values = custom_palette)

gene_signature <- as.data.frame(cbind(bk.meta$Run,as.numeric(bk.meta$AC_UCell),
                                      as.numeric(bk.meta$Neu_UCell),
                                      as.numeric(bk.meta$OL_UCell),
                                      as.numeric(bk.meta$OPC_UCell),
                                      as.numeric(bk.meta$BC_UCell),
                                      as.numeric(bk.meta$TC_UCell),
                                      as.numeric(bk.meta$NK_UCell),
                                      as.numeric(bk.meta$PVM_UCell),
                                      as.numeric(bk.meta$Mono_UCell),
                                      as.numeric(bk.meta$cDC_UCell),
                                      as.numeric(bk.meta$MG_UCell),
                                      as.numeric(bk.meta$pDC_UCell),
                                      as.numeric(bk.meta$EC_UCell),
                                      as.numeric(bk.meta$FBMC_UCell),
                                      as.numeric(bk.meta$aFB2_UCell),
                                      as.numeric(bk.meta$FB_UCell),
                                      as.numeric(bk.meta$aFB1_UCell),
                                      as.numeric(bk.meta$aMac_UCell),
                                      as.numeric(bk.meta$SMC_UCell)
))
seed_run <- bk.meta$Run[which(bk.meta$cell_color2 == 0)]
colnames(gene_signature) <- c('Run','AC','Neu','OL','OPC',
                              'BC','TC','NK','PVM','Mono','cDC','MG',
                              'pDC','EC','FBMC','aFB2','FB',
                              'aFB1','aMac','SMC'
)
gene_signature_values <- gene_signature[,2:ncol(gene_signature)]
seed_idx <- which(gene_signature$Run == seed_run)
for(i in 1:ncol(gene_signature_values)){
  setpoint <- as.numeric(gene_signature_values[seed_idx,i])
  gene_signature_values[,i] <- as.numeric(gene_signature_values[,i]) - setpoint
}
gene_signature <- cbind(gene_signature$Run,gene_signature_values)
colnames(gene_signature)[1] <- 'Run'
gene_signature %<>% 
  as_tibble() %>% 
  pivot_longer(cols = c('AC','Neu','OL','OPC',
                        'BC','TC','NK','PVM','Mono','cDC','MG',
                        'pDC','EC','FBMC','aFB2','FB',
                        'aFB1','aMac','SMC'
  ),names_to = 'celltype')
bk.meta %>% 
  left_join(gene_signature,by = c('Run')) %>%
  ggplot(aes(cell_color2,(as.numeric(value)))) +
  geom_point(aes(color = celltype)) + 
  geom_smooth(aes(color = celltype),method='loess',se = F) + 
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_color_manual(values = custom_palette)

# ANOVA analysis with plots
bk.seurat <-  readRDS('Final_bulk_dx_pseudotime_obj.RDS')
custom_palette_tissue.status <- c('grey','darkorange','red')
names(custom_palette_tissue.status) <- c('Intracranial cortical artery','Unruptured intracranial aneurysm','Ruptured intracranial aneurysm')
bk.seurat@meta.data %>% 
  ggplot(aes(x = factor(tissue.status,levels = c('Intracranial cortical artery','Unruptured intracranial aneurysm','Ruptured intracranial aneurysm')),
             y = cell_color2,
             fill = tissue.status,
  )) +
  geom_boxplot(outliers = F) +
  theme_classic() +
  geom_jitter(width = 0.1,alpha = 0.25,size = 3) + 
  scale_fill_manual(values =custom_palette_tissue.status)
aov.results <- aov(cell_color2 ~ tissue.status,data = bk.seurat@meta.data) #stats
summary(aov.results)
tukey.two.way <- TukeyHSD(aov.results) #stats
tukey.two.way



