obj <- readRDS(file = "Final_ax_object.RDS")

#Deconvolutions on pseudobulks
CTRL <- obj
pseudo.CTRL <- AggregateExpression(CTRL, assays = "SCT", return.seurat = T, group.by = c("Sample"))
pseudo.CTRL.matrix <- as.matrix(t(pseudo.CTRL@assays[["SCT"]]@layers[["counts"]]))
rownames(pseudo.CTRL.matrix) <- colnames(pseudo.CTRL)
colnames(pseudo.CTRL.matrix) <- rownames(pseudo.CTRL)
bk.dat <- pseudo.CTRL.matrix

sc.dat <- obj
sc.dat <- t(as.data.frame(sc.dat@assays[["SCT"]]@counts))
cell.state.labels <- as.character(obj$filtered_annotations)
cell.type.labels <- cell.state.labels
cell.type.labels[which(cell.type.labels %in% c('PVM','Mono','aMac','cDC','MG'))] <- 'Myeloid'

table(cbind.data.frame(cell.state.labels, cell.type.labels))

plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell type correlation",
              cexRow=0.7, cexCol=0.7,
              margins=c(4,4))

sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=cell.type.labels,
  species="hs",
  return.raw=TRUE
)

head(sc.stat)  

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=sc.dat, 
  cell.type.labels=cell.type.labels,
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

theta <- get.fraction (bp=bp.res,
                       which.theta="first",
                       state.or.type="state")
head(theta)
theta.cv <- bp.res@posterior.theta_f@theta.cv
head(theta.cv)

theta <- as.data.frame(theta)
theta$Run <- rownames(theta)
theta <- pivot_longer(theta,cols = -c("Run"))

props <- getTransformedProps(clusters = CTRL$filtered_annotations, 
                             sample = CTRL$Sample)
index <- rowSums(props[["Counts"]])>0
props[["Proportions"]]  <- props[["Proportions"]][index,]
real.full <- props$Proportions
real.full %<>% 
  as_tibble() 
colnames(real.full) <- c('name','Run','real.value')

converter <- cbind(setdiff(real.full$Run,theta$Run),setdiff(theta$Run,real.full$Run))

colnames(theta)[3] <- 'simulated.value'
for(i in unique(theta$Run)){
  print(i)
  if(i %in% converter[,2]){
    index <- which(converter[,2] == i)
    y <- converter[index,1]
  }else{
    y <- i
  }
  simulated <- theta %>% 
    filter(Run == i) %>% 
    mutate(Run = y)
  real <- real.full %>% 
    filter(Run == y) 
  merged <- inner_join(simulated,real,by = c('name','Run'))
  R <- cor(merged$simulated.value,merged$real.value)
  if(i == unique(theta$Run)[1]){
    merged.full <- merged
    R.full <- R
  }else{
    merged.full <- rbind(merged.full,merged)
    R.full <- c(R.full,R)
  }
}
merged.full.reserve <- merged.full

#Deconvolutions on bulk-RNAseq samples
raw_counts <- read.csv(file = "Final_bulkcontrols_sleuth.csv",row.names = 1)
obj <- readRDS(file = "Final_ax_object.RDS")
CTRL <- obj
bk.dat <- raw_counts


sc.dat <- obj
sc.dat <- t(as.data.frame(sc.dat@assays[["SCT"]]@counts))
cell.state.labels <- as.character(obj$filtered_annotations)
cell.type.labels <- cell.state.labels
cell.type.labels[which(cell.type.labels %in% c('PVM','Mono','aMac','cDC','MG'))] <- 'Myeloid'

table(cbind.data.frame(cell.state.labels, cell.type.labels))

plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell type correlation",
              cexRow=0.7, cexCol=0.7,
              margins=c(4,4))

sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=cell.type.labels,
  species="hs",
  return.raw=TRUE 
)

head(sc.stat)  

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=sc.dat,  
  cell.type.labels=cell.type.labels,
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

theta <- get.fraction (bp=bp.res,
                       which.theta="first",
                       state.or.type="state")
head(theta)
theta.cv <- bp.res@posterior.theta_f@theta.cv
head(theta.cv)

theta <- as.data.frame(theta)
theta$Run <- rownames(theta)
theta <- pivot_longer(theta,cols = -c("Run"))

Idents(CTRL) <- CTRL$Sample
CTRL_subset <- CTRL %>% 
  subset(idents = c('ctrl085_L','ctrl086_1','ctrl086_2','ctrl099_21','ctrl099_22'))
CTRL_subset$Sample_collapsed <- CTRL_subset$Sample
CTRL_subset$Sample_collapsed[which(CTRL_subset$Sample_collapsed %in% c('ctrl086_1','ctrl086_2'))] <- 'ctr86'
CTRL_subset$Sample_collapsed[which(CTRL_subset$Sample_collapsed %in% c('ctrl085_L'))] <- 'ctr85'
CTRL_subset$Sample_collapsed[which(CTRL_subset$Sample_collapsed %in% c('ctrl099_21'))] <- 'ctr9921'
CTRL_subset$Sample_collapsed[which(CTRL_subset$Sample_collapsed %in% c('ctrl099_22'))] <- 'ctr9922'
props <- getTransformedProps(clusters = CTRL_subset$filtered_annotations, 
                             sample = CTRL_subset$Sample_collapsed)
index <- rowSums(props[["Counts"]])>0
props[["Proportions"]]  <- props[["Proportions"]][index,]
real.full <- props$Proportions
real.full %<>% 
  as_tibble() 
colnames(real.full) <- c('name','Run','real.value')

converter <- cbind(unique(real.full$Run),unique(theta$Run))

colnames(theta)[3] <- 'simulated.value'
for(i in unique(theta$Run)){
  print(i)
  if(i %in% converter[,2]){
    index <- which(converter[,2] == i)
    y <- converter[index,1]
  }else{
    y <- i
  }
  simulated <- theta %>% 
    filter(Run == i) %>% 
    mutate(Run = y)
  real <- real.full %>% 
    filter(Run == y) 
  merged <- inner_join(simulated,real,by = c('name','Run'))
  R <- cor(merged$simulated.value,merged$real.value)
  if(i == unique(theta$Run)[1]){
    merged.full <- merged
    R.full <- R
  }else{
    merged.full <- rbind(merged.full,merged)
    R.full <- c(R.full,R)
  }
}
colnames(merged.full)[3] <- 'bulk.value'
colnames(merged.full)[4] <- 'singlecell.value'


