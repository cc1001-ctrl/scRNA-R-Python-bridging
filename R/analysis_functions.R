# Basic single-cell RNA-seq analysis functions
# Author: CC1001-CTRL
# Date: 2024

#' Load and preprocess single-cell data
#'
#' @param data_path Path to the data file
#' @param project_name Name for the Seurat object
#' @return A preprocessed Seurat object
#' @export
load_and_preprocess <- function(data_path, project_name = "scRNA_analysis") {
  library(Seurat)
  library(dplyr)
  
  # Load data (adjust based on your data format)
  # This is a template - modify for your specific data format
  if (file.exists(data_path)) {
    data <- Read10X(data_path)
    seurat_obj <- CreateSeuratObject(counts = data, project = project_name)
  } else {
    # Use example data if file doesn't exist
    data("pbmc_small")
    seurat_obj <- pbmc_small
    message("Using example pbmc_small dataset")
  }
  
  # Basic preprocessing
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells and features
  seurat_obj <- subset(seurat_obj, 
                      subset = nFeature_RNA > 200 & 
                               nFeature_RNA < 5000 & 
                               percent.mt < 20)
  
  return(seurat_obj)
}

#' Perform standard normalization and scaling
#'
#' @param seurat_obj A Seurat object
#' @param nfeatures Number of variable features to find
#' @return A normalized and scaled Seurat object
#' @export
normalize_and_scale <- function(seurat_obj, nfeatures = 2000) {
  library(Seurat)
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                    nfeatures = nfeatures)
  
  # Scale data
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  return(seurat_obj)
}

#' Perform dimensional reduction
#'
#' @param seurat_obj A Seurat object
#' @param npcs Number of principal components to compute
#' @return A Seurat object with PCA and UMAP
#' @export
run_dimensional_reduction <- function(seurat_obj, npcs = 50) {
  library(Seurat)
  
  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj),
                      npcs = npcs)
  
  # Find neighbors and clusters
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  
  return(seurat_obj)
}

#' Create quality control plots
#'
#' @param seurat_obj A Seurat object
#' @return A list of ggplot objects
#' @export
create_qc_plots <- function(seurat_obj) {
  library(ggplot2)
  library(Seurat)
  
  plots <- list()
  
  # Violin plot for QC metrics
  plots$violin <- VlnPlot(seurat_obj, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                         ncol = 3)
  
  # Feature scatter plots
  plots$feature_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", 
                                          feature2 = "percent.mt")
  plots$feature_scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", 
                                          feature2 = "nFeature_RNA")
  
  return(plots)
}