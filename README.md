# scRNA-R-Python-bridging

A minimal R single-cell template for single-cell RNA sequencing analysis with Seurat, SpatialExperiment, and Bioconductor packages.

## Quick Start

1. **Clone this repository:**
   ```bash
   git clone https://github.com/cc1001-ctrl/scRNA-R-Python-bridging.git
   cd scRNA-R-Python-bridging
   ```

2. **Install R dependencies using renv:**
   ```r
   # Install renv if not already installed
   if (!requireNamespace("renv", quietly = TRUE)) {
     install.packages("renv")
   }
   
   # Restore the package environment
   renv::restore()
   ```

3. **Load required libraries:**
   ```r
   library(Seurat)
   library(SpatialExperiment)
   library(ggplot2)
   library(dplyr)
   ```

## Project Structure

```
├── R/                  # R scripts and functions
├── data/               # Input data files
├── outputs/            # Analysis outputs and results
├── vignettes/          # Analysis notebooks and tutorials
├── renv.lock           # R package dependencies
└── README.md           # This file
```

## Dependencies

This template includes the following key packages:

- **Seurat** (v4.3.0): Comprehensive single-cell RNA-seq analysis
- **SpatialExperiment** (v1.10.0): Spatial transcriptomics data handling
- **ggplot2** (v3.4.2): Data visualization
- **Bioconductor** packages: Core infrastructure for genomics data
- **dplyr**: Data manipulation
- **Matrix**: Sparse matrix operations

## renv Package Management

This project uses `renv` for reproducible package management. The `renv.lock` file contains all package versions and dependencies.

### Initial Setup
```r
# Install renv
install.packages("renv")

# Initialize renv in your project (if not already done)
renv::init()

# Restore packages from lockfile
renv::restore()
```

### Managing Dependencies
```r
# Install new packages
install.packages("package_name")
renv::snapshot()  # Update lockfile

# Update all packages
renv::update()

# Check package status
renv::status()
```

## Example Usage

```r
# Load libraries
library(Seurat)
library(ggplot2)

# Create a simple Seurat object example
# (Replace with your actual data loading)
data("pbmc_small")
pbmc <- pbmc_small

# Basic analysis workflow
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualization
DimPlot(pbmc, reduction = "pca")
```

## Getting Help

- [Seurat Documentation](https://satijalab.org/seurat/)
- [SpatialExperiment Documentation](https://bioconductor.org/packages/SpatialExperiment/)
- [Bioconductor Documentation](https://bioconductor.org/)

## Citation

If you use this template in your research, please cite:

```bibtex
@software{scRNA_R_Python_bridging,
  title = {scRNA-R-Python-bridging: Minimal R Single-Cell Template},
  author = {{CC1001-CTRL}},
  year = {2024},
  url = {https://github.com/cc1001-ctrl/scRNA-R-Python-bridging},
  note = {R template for single-cell RNA sequencing analysis}
}
```

### Package Citations

- **Seurat**: Hao et al. "Integrated analysis of multimodal single-cell data." Cell (2021)
- **SpatialExperiment**: Righelli et al. "SpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor." Bioinformatics (2022)
- **ggplot2**: Wickham, H. "ggplot2: Elegant Graphics for Data Analysis." Springer-Verlag New York (2016)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
