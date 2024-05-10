# Propeller on CTRL v AX
obj <- readRDS(file='Final_ax_object.RDS')
DF <- obj
DF <- subset(DF,idents = names(table(DF$filtered_annotations))[-which(names(table(DF$filtered_annotations)) %in% c("Neu","AC","OPC","OL"))])

Ordered <- factor(DF$Sample, 
                  levels=c('a1217','a1211','Ax-12-13','Ax-3623','Ax2623_TB039','TB045','Ax1116','Ax4323','ctrl086_2','ctrl086_1','ctrl099_21',
                           'ctrl085_L','UCSF-2023-003','UCSF-2022-006','UCSF-2021-009','ctrl099_22'))

p1 <- plotCellTypeProps(clusters = DF$filtered_annotations, sample = Ordered) + theme(axis.text.x = element_text(angle = 45))+ ggtitle("Ax Cell Type Proportions") + 
  theme(plot.title = element_text(size = 18, hjust = 0)) +
  scale_fill_manual(values=custom_palette)
p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90))
props <- getTransformedProps(clusters = DF$filtered_annotations, 
                             sample = DF$Sample)
index <- rowSums(props[["Counts"]])>0
props[["Counts"]] <- props[["Counts"]][index,]
props[["TransformedProps"]]  <- props[["TransformedProps"]][index,]
props[["Proportions"]]  <- props[["Proportions"]][index,]
grp <- c(rep("Unruptured_Unstable",7),rep("CTRL",5),'Unruptured_Unstable',rep("CTRL",3))
designAS <- model.matrix(~0+grp)
colnames(designAS) <- c("CTRL","Unruptured_Unstable")
mycontr <- makeContrasts(Unruptured_Unstable-CTRL, levels=designAS)
results <- propeller.ttest(prop.list = props,design = designAS, contrasts = mycontr,
                           robust=TRUE,trend=FALSE,sort=TRUE)
results$significance <- "black"
results$significance[which(results$FDR<0.05)] <- "red"


