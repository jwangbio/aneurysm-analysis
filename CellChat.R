rm(list=ls())
.rs.restartR()
set.seed(67)
gc()


library(CellChat)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 891289600)

obj <- readRDS('Final_ax_object.RDS')
CTRL <- subset(obj,subset = Condition == 'control')
data.input <- CTRL[["SCT"]]@data # normalized data matrix
labels <- Idents(CTRL)
meta <- data.frame(labels = labels, row.names = names(labels))
CTRL$samples <- factor(CTRL$Sample)

cellchat <- createCellChat(object = CTRL, group.by = "ident", assay = "SCT")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 120) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = T,population.size = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,'Final_cellchat_CTRL.RDS')

obj <- readRDS('Final_ax_object.RDS')
AX <- subset(obj,subset = Condition == 'aneurysm')
data.input <- AX[["SCT"]]@data # normalized data matrix
labels <- Idents(AX)
meta <- data.frame(labels = labels, row.names = names(labels))
AX$samples <- factor(AX$Sample)

cellchat <- createCellChat(object = AX, group.by = "ident", assay = "SCT")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 120) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = T,population.size = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,'Final_cellchat_AX.RDS')

# Comparing the two
cellchat_CTRL <- readRDS('Final_cellchat_CTRL.RDS')
cellchat_AX <- readRDS('Final_cellchat_AX.RDS')
object.list <- list(CTRL = cellchat_CTRL,AX = cellchat_AX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'aFB1')
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'aMac',)
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'SMC',)
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'aFB1',measure = 'weight')
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'aMac',measure = 'weight')
netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = 'SMC',measure = 'weight')
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


object.list <- list(CTRL = cellchat_CTRL,AX = cellchat_AX)
group.cellType <- c(rep("aFB1", 2), rep("aMac", 2), rep("cDC", 2),rep("Mono", 2),rep("PVM", 2),rep("MG", 2))
group.cellType <- factor(group.cellType, levels = c("aFB1", "aMac",'cDC','Mono','PVM','MG'))
custom_palette_FB_Myeloid <- custom_palette[c("aFB1", "aMac",'cDC','Mono','PVM','MG')]

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat2 <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "count.merged", label.edge = F,color.use = custom_palette_FB_Myeloid)
netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "weight.merged", label.edge = F,color.use = custom_palette_FB_Myeloid)

gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = custom_palette,cluster.rows = T,cluster.cols = T)
#gg1 +
#gg2

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# Custom differential in/out by celltype
pg <- ggplot_build(gg[[1]]) 
pg_CTRL <- pg$data[[2]][,c('label','x','y')]
pg <- ggplot_build(gg[[2]])
pg_AX <- pg$data[[2]][,c('label','x','y')]
Count <- num.link[,2]
colnames(pg_CTRL)[2:3] <- c('out_ctrl','in_ctrl')
colnames(pg_AX)[2:3] <- c('out_ax','in_ax')
pg <- left_join(pg_CTRL,pg_AX,by = 'label')
pg$Count <- Count
pg %<>% 
  mutate(out_diff = out_ax - out_ctrl,
         in_diff = in_ax - in_ctrl
  )
ggplot(pg,aes(out_diff,in_diff)) +
  geom_point(aes(size = Count,color = label),) +
  geom_text_repel(aes(label = label)) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic() + 
  scale_color_manual(values = custom_palette)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "aMac")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "aFB1")
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,color.use = c('darkgrey','#ff001b'))
#gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 #+ gg2

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height =16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
#for (i in 1:length(object.list)) {
for (i in 1) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("GALECTIN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
#for (i in 1:length(object.list)) {
for (i in 2) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


pathways.show <- c("CXCL") 
 vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = 'MIF', layout = "chord")
netVisual_chord_gene(cellchat, sources.use = 'aFB1', targets.use = c('aMac','cDC'), lab.cex = 0.5,legend.pos.y = 30,slot.name = 'netP')
netVisual_heatmap(cellchat, signaling = pathways.show)

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram
group.cellType <- levels(object.list[[1]]@idents)
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
par(mfrow = c(1, 1), xpd=TRUE)
for (i in 2) {
  netVisual_chord_gene(object.list[[i]], sources.use = 'aFB1',
                       targets.use = c('aMac','cDC','Mono','PVM','MG'),
                       lab.cex = 0.5, 
                       color.use = custom_palette_FB_Myeloid
  )
}
par(mfrow = c(1, 1), xpd=TRUE)
for (i in 2) {
  netVisual_chord_gene(object.list[[i]], targets.use = 'aFB1',
                       sources.use = c('aMac','cDC','Mono','PVM','MG'),
                       lab.cex = 0.5, 
                       color.use = custom_palette_FB_Myeloid
  )
}







