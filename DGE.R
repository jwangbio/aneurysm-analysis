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
ibrary(ggrepel)
library(apeglm)

#Looped DESeq2
obj <- readRDS('Final_ax_object.RDS')
pseudobulk <- AggregateExpression(obj, assays = "SCT", group.by = c('Sample',"filtered_annotations",'ax_status','Preparation'), return.seurat = TRUE)
all_res_df <- data.frame()
cell_types <- c('FB','SMC','Art','FBMC')

for (cell_type in cell_types) {
  message(cell_type)
  
  # Subset to current cell type
  cell_subset <- subset(pseudobulk, subset = filtered_annotations == cell_type)
  counts <- as.matrix(cell_subset@assays$SCT$counts)
  coldata <- cell_subset@meta.data
  
  # Ensure ax_status is a factor with CTRL as reference
  coldata$ax_status <- factor(coldata$ax_status, levels = c("CTRL","Unruptured"))
  
  # Preparation as factor
  coldata$Preparation <- factor(coldata$Preparation)
  
  # Keep only covariates with >1 level
  covariates <- c("ax_status", "Preparation")
  valid_covariates <- covariates[sapply(covariates, function(cov) length(unique(coldata[[cov]])) > 1)]
  
  # Check that ax status is >1 level 
  if (!"ax_status" %in% valid_covariates) {
    message("Skipping: ax_status has only one level for this cell type.")
    next
  }
  
  # Design formula
  design_formula <- as.formula(paste("~", paste(valid_covariates, collapse = " + ")))
  
  # DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_formula)
  
  # Filter low-expressed genes
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run DESeq with Wald test and shrinkage
  dds <- DESeq(dds, test = "Wald")
  res <- results(dds, contrast = c("ax_status", "Unruptured", "CTRL"))
  res <- lfcShrink(dds, coef = "ax_status_Unruptured_vs_CTRL", type = "apeglm")
  
  # Prepare output
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- cell_type
  
  # Significance thresholds
  log2FC_threshold <- 1
  padj_threshold <- 0.05
  res_df$significant <- ifelse(!is.na(res_df$padj) &
                                 res_df$padj < padj_threshold &
                                 abs(res_df$log2FoldChange) > log2FC_threshold,
                               "Significant", "Not significant")
  
  # Accumulate
  all_res_df <- rbind(all_res_df, res_df)
}
