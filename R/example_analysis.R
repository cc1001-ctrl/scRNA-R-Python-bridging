# Example single-cell RNA-seq analysis workflow
# This script demonstrates a complete analysis pipeline

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source custom functions
source("R/analysis_functions.R")

# Set seed for reproducibility
set.seed(42)

# 1. Load and preprocess data
message("Loading and preprocessing data...")
seurat_obj <- load_and_preprocess("data/")  # Will use example data if path doesn't exist

# 2. Normalize and scale data
message("Normalizing and scaling data...")
seurat_obj <- normalize_and_scale(seurat_obj, nfeatures = 2000)

# 3. Run dimensional reduction
message("Running dimensional reduction...")
seurat_obj <- run_dimensional_reduction(seurat_obj, npcs = 50)

# 4. Create visualizations
message("Creating visualizations...")

# Quality control plots
qc_plots <- create_qc_plots(seurat_obj)

# Save QC plots
ggsave("outputs/qc_violin_plot.png", qc_plots$violin, width = 12, height = 6)
ggsave("outputs/qc_scatter_plot1.png", qc_plots$feature_scatter1, width = 8, height = 6)
ggsave("outputs/qc_scatter_plot2.png", qc_plots$feature_scatter2, width = 8, height = 6)

# PCA plot
pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters")
ggsave("outputs/pca_plot.png", pca_plot, width = 8, height = 6)

# UMAP plot
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", 
                     label = TRUE, pt.size = 0.5) + 
             ggtitle("UMAP Clustering")
ggsave("outputs/umap_plot.png", umap_plot, width = 8, height = 6)

# Variable features plot
var_features_plot <- VariableFeaturePlot(seurat_obj)
top10_features <- head(VariableFeatures(seurat_obj), 10)
var_features_plot <- LabelPoints(plot = var_features_plot, points = top10_features, repel = TRUE)
ggsave("outputs/variable_features_plot.png", var_features_plot, width = 10, height = 6)

# 5. Find marker genes for clusters
message("Finding marker genes...")
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                             logfc.threshold = 0.25)

# Save marker genes
write.csv(all_markers, "outputs/cluster_markers.csv", row.names = FALSE)

# Plot top markers
top_markers <- all_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)

if(nrow(top_markers) > 0) {
  marker_plot <- FeaturePlot(seurat_obj, features = head(top_markers$gene, 6), 
                            ncol = 3)
  ggsave("outputs/top_markers_plot.png", marker_plot, width = 15, height = 10)
}

# 6. Create summary statistics
message("Creating summary statistics...")
summary_stats <- data.frame(
  n_cells = ncol(seurat_obj),
  n_genes = nrow(seurat_obj),
  n_clusters = length(unique(Idents(seurat_obj))),
  median_genes_per_cell = median(seurat_obj$nFeature_RNA),
  median_UMI_per_cell = median(seurat_obj$nCount_RNA),
  median_mito_percent = median(seurat_obj$percent.mt)
)

write.csv(summary_stats, "outputs/analysis_summary.csv", row.names = FALSE)

# Save the final Seurat object
saveRDS(seurat_obj, "outputs/processed_seurat_object.rds")

message("Analysis complete! Check the outputs/ directory for results.")
print(summary_stats)