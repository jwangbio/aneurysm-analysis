# Pseudotime on Annoted FB
custom_palette_FB <- unique(c("#e5e1e3",#aFB2
                              "#ff001b",#aFB1
                              "#fbe29e",#FB1
                              '#b9ceff',#FB2
                              "#0085ff"#FB3
))
names(custom_palette_FB) <- c('aFB2','aFB1','FB1','FB2','FB3')
subset <- readRDS("Final_FB_object.RDS")
cds <- SeuratWrappers::as.cell_data_set(subset)
cds <- cluster_cells(cds,partition_qval = 1,random_seed = 67)
plot_cells(cds, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
           color_cells_by = 'fb_subclusters',show_trajectory_graph = F,label_cell_groups = F,rasterize = F) +
  scale_color_manual(values = custom_palette_FB)

max <- which.max(unlist(FetchData(subset, "AVP")))
max <- colnames(subset)[max]
cds <- order_cells(cds, root_cells = max)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE,show_trajectory_graph = F,rasterize = F) + 
  scale_color_viridis()  

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=120)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value == 0 ))
ciliated_cds_pr_test_res %>% 
  mutate(X = rownames(.)) %>% 
  slice_max(n = 15,order_by = morans_test_statistic) %>% 
  ggplot(mapping = aes(x = morans_test_statistic,y = reorder(X,morans_test_statistic)))+
  geom_bar(stat="identity",) +
  ylab("Top 15 Genes")

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

cds <- preprocess_cds(cds)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], 
                                    resolution=4e-2,
                                    random_seed = 67,
)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$fb_subclusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
all_modules <- gene_module_df %>% as.data.frame()
module_explore <- gene_module_df %>% as.data.frame() %>% filter(module %in% c(8,16)) #aFB1  at 4e-2       SELECTED RNG = 67   
module_explore <- gene_module_df %>% as.data.frame() %>% filter(module %in% c(22,24)) #aFB2 at 4e-2       SELECTED RNG = 67
module_explore_filtered <- module_explore %>% filter(!(str_detect(id, "^ENSG"))) %>% filter(!(str_detect(id, "^MT-"))) %>% filter(!(str_detect(id, "^RPS"))) %>% filter(!(str_detect(id, "^RPL")))

module_signature <- list(
  module_genes = module_explore_filtered$id
)
aneurysm_module <-AddModuleScore_UCell(subset,features=module_signature, name= NULL)
aneurysm_module$pseudotime <- pseudotime(cds)
Idents(aneurysm_module) <- 'NA'
FeatureScatter(aneurysm_module, feature1 = "pseudotime", feature2 = names(module_signature), pt.size = 1, jitter = TRUE, seed = 67, 
               cols='grey'
) + 
  geom_smooth(method = "loess", se = F, span = 0.95, linewidth=2,col = 'black')

