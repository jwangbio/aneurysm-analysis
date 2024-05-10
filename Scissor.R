obj <- readRDS(file='Final_ax_object.RDS')
DF <- obj
DF[['cleaned']] <- NULL
DF <- RenameAssays(DF,assay.name = 'SCT',new.assay.name = 'RNA')
DF@graphs[["RNA_nn"]] <- DF@graphs[["SCT_nn"]]
DF@graphs[["RNA_snn"]] <- DF@graphs[["SCT_snn"]]
DF@graphs[["SCT_snn"]] <- NULL
DF@graphs[["SCT_nn"]] <- NULL

bk.dat <- read.csv("Final_bulk_counts_matrix.csv",row.names = NULL)
duplicates <- bk.dat$row.names[duplicated(bk.dat$row.names)]
bk.dat <- bk.dat[-which(bk.dat$row.names %in% duplicates),]
rownames(bk.dat) <- bk.dat$row.names
bk.dat <- bk.dat[,-1]
bk.meta <- read.csv("Final_geo_aneurysm_meta_data.csv")
bk.meta %<>% 
  as_tibble() %>%
  filter(Run %in% colnames(bk.dat)) %>% 
  filter(rupture.status %in% c('Ruptured','Unruptured'))
bk.dat <- bk.dat[,bk.meta$Run]
all(colnames(bk.dat) == bk.meta$Run)
phenotype <- bk.meta$rupture.status
phenotype[which(phenotype == 'Ruptured')] <- 1
phenotype[which(phenotype == 'Unruptured')] <- 0
phenotype <- as.numeric(phenotype)
tags <- unique(bk.meta$rupture.status)
infos1 <- Scissor2(bk.dat, DF, phenotype, alpha = NULL, cutoff = 0.08,tag = tags,
                   family = "binomial")

# Post Scissor analysis
obj$Scissor3 <- NA
obj$Scissor3[which(obj$Cell_ID %in% infos1$Scissor_pos)] <- 'Rupture'
obj$Scissor3[which(obj$Cell_ID %in% infos1$Scissor_neg)] <- 'Unruptured'
obj$Scissor_alpha <- 0.05
obj$Scissor_alpha[!is.na(obj$Scissor3)] <- 1
Scissor_colors <- c('Rupture','Unruptured','NA')
Scissor_colors <- c('red','darkorange','lightgrey')

#Stats on Scissor Celltype - Hypergeometric Test
meta <- obj@meta.data
meta %<>% 
  as_tibble() %>% 
  mutate(ruptured_binary = case_when(Scissor3 == 'Rupture' ~ 1,
                                     is.na(Scissor3) | Scissor3 != 'Rupture' ~ 0),
         unruptured_binary = case_when(Scissor3 == 'Unruptured' ~ 1,
                                       is.na(Scissor3) | Scissor3 != 'Unruptured'  ~ 0)
  )

pval <- c() #ruptured 
FC <- c()
for(i in unique(meta$filtered_annotations)){
  print(i)
  sample_n <- as.numeric(table(meta$filtered_annotations)[which(names(table(meta$filtered_annotations)) == i)])
  sample_hits <- as.numeric(table(meta$filtered_annotations,meta$ruptured_binary)[i,2])
  pop_n <- nrow(meta)
  pop_hits <- as.numeric(table(meta$ruptured_binary)[2])
  expected_ratio <- pop_hits/pop_n
  sample_ratio <- sample_hits/sample_n
  pval <- c(pval,phyper(sample_hits,pop_hits,pop_n - pop_hits,sample_n,lower.tail = FALSE))
  FC <- c(FC,sample_ratio / expected_ratio)
}
names(pval) <- unique(meta$filtered_annotations)
fdr <- p.adjust(pval,method = 'fdr')
names(FC) <- unique(meta$filtered_annotations)
FC[which(pval < 0.05)]

pval <- c() #unruptured
FC <- c()
for(i in unique(meta$filtered_annotations)){
  print(i)
  sample_n <- as.numeric(table(meta$filtered_annotations)[which(names(table(meta$filtered_annotations)) == i)])
  sample_hits <- as.numeric(table(meta$filtered_annotations,meta$unruptured_binary)[i,2])
  pop_n <- nrow(meta)
  pop_hits <- as.numeric(table(meta$unruptured_binary)[2])
  expected_ratio <- pop_hits/pop_n
  sample_ratio <- sample_hits/sample_n
  pval <- c(pval,phyper(sample_hits,pop_hits,pop_n - pop_hits,sample_n,lower.tail = FALSE))
  FC <- c(FC,sample_ratio / expected_ratio)
}
names(pval) <- unique(meta$filtered_annotations)
fdr <- p.adjust(pval,method = 'fdr')
names(FC) <- unique(meta$filtered_annotations)
FC[which(fdr < 0.05)]


