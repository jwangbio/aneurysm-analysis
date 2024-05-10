rm(list=ls())
.rs.restartR()
set.seed(67)
gc()

library(Seurat)
library(UCell)
library(purrr)

CTRL <- readRDS('Final_CTRL_spatial.rds')
DefaultAssay(CTRL) <- 'RNA'
AX <- readRDS('Final_AX_spatial.rds')
DefaultAssay(AX) <- 'RNA'

SpatialFeaturePlot(CTRL, features = 'POSTN',max.cutoff = 5,pt.size.factor = 3,alpha = c(0.05,1))
SpatialFeaturePlot(AX, features = 'POSTN',max.cutoff = 5,pt.size.factor = 3,alpha = c(0.05,1))

# Cell marker specific signatures
DefaultAssay(AX) <- 'RNA'
DefaultAssay(CTRL) <- 'RNA'
gene_list <- read.csv('Final_annotated_markers.csv')
gene_list %<>% 
  group_by(cluster) %>%
  slice_max(n = 20,order_by = weighted_specificity) %>% 
  ungroup() %>% 
  select(cluster,gene)
spatial.intersect <- intersect(rownames(CTRL),rownames(AX))
intersect <- intersect(gene_list$gene,spatial.intersect)
gene_list %<>% 
  filter(gene %in% intersect)
Tgenes_all <- c()
for(i in unique(gene_list$cluster)){
  print(i)
  Tgenes <- gene_list %>% 
    filter(cluster == i) %>% 
    pull(gene)
  Tgenes_all <- c(Tgenes_all,list(Tgenes))
}
names(Tgenes_all) <- unique(gene_list$cluster)
CTRL <- AddModuleScore_UCell(CTRL, features=Tgenes_all)
AX <- AddModuleScore_UCell(AX, features=Tgenes_all) 

AX$condition <- 'AX'
CTRL$condition <- 'CTRL'
merged <- merge(CTRL,AX)

# FB Activation
custom_palette_ax_status <- unique(c(
  "#e5e1e3",#CTRL
  "#ff001b"#Unruptured
))
names(custom_palette_ax_status) <- c('CTRL','AX')
fb_activation <- read.table('Final_HCM_activatedFB_gene_list.txt',header = T)
fb_activation <- as.data.frame(fb_activation[1:20,])
colnames(fb_activation) <- 'gene'
spatial.intersect <- intersect(fb_activation$gene,rownames(merged))
intersect <- intersect(fb_activation$gene,spatial.intersect)
gene_sig <- list(fb_activation = intersect)
merged <- AddModuleScore_UCell(merged,features = gene_sig,ncores = 120)
merged2 <- subset(merged,subset = condition == 'AX' & nCount_RNA > 1000 | condition == 'CTRL' & nCount_RNA > 1000)
merged2@meta.data %>% 
  pivot_longer(cols = c('fb_activation_UCell'),names_to = 'gene_set') %>% 
  ggplot(aes(x=factor(condition,levels = c('CTRL','AX')),y = value,fill = factor(condition,levels = c('CTRL','AX')))) +
  geom_boxplot() +
  facet_grid(cols = vars(gene_set)) +
  theme_classic() +
  scale_fill_manual(values = custom_palette_ax_status)
SpatialFeaturePlot(merged, features = 'fb_activation_UCell',max.cutoff = 'q99',keep.scale = 'all',pt.size.factor = 2,alpha = c(0.05,1))

# VlnPlots
custom_palette_ax_status <- unique(c(
  "#e5e1e3",#CTRL
  "#ff001b"#Unruptured
))
names(custom_palette_ax_status) <- c('CTRL','AX')
merged$condition <- factor(merged$condition,levels = c('CTRL','AX'))
VlnPlot(merged, features = c("log10_nCount_RNA"), pt.size = 0.0,group.by = 'condition') + NoLegend() +
  scale_fill_manual(values = custom_palette_ax_status)
t_test(merged@meta.data,formula = log10_nCount_RNA ~ condition)
VlnPlot(merged, features = c("log10_nFeature_RNA"), pt.size = 0.0,group.by = 'condition') + NoLegend() +
  scale_fill_manual(values = custom_palette_ax_status)
t_test(merged@meta.data,formula = log10_nFeature_RNA ~ condition)

subset <- subset(merged,subset = condition == 'AX' & nCount_RNA > 1000 | condition == 'CTRL' & nCount_RNA > 1000)
SpatialFeaturePlot(subset, features = c('aMac_UCell'
),keep.scale = 'all',pt.size.factor = 2,alpha = c(0.05,1))
subset$condition <- factor(subset$condition,levels = c('CTRL','AX'))
subset@meta.data %>% 
  pivot_longer(cols = c('PVM_UCell',
  ),names_to = 'gene_set') %>% 
  ggplot(aes(x=condition,y = value,fill = condition)) +
  geom_boxplot(outliers = F) +
  facet_grid(cols = vars(gene_set)) +
  scale_fill_manual(values = custom_palette_ax_status) + 
  theme_classic()
VlnPlot(merged, features = c('aMac_UCell','aFB1_UCell','SMC_UCell',
                             'BC_UCell','FB_UCell','Mono_UCell',
                             'aFB2_UCell','FBMC_UCell','NK_UCell',
                             'cDC_UCell','TC_UCell','EC_UCell'
), pt.size = 0.0,group.by = 'condition') + NoLegend()




