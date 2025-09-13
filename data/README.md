# Data Directory

This directory is intended for input data files for single-cell RNA sequencing analysis.

## Recommended Data Formats

### 10X Genomics Output
- `barcodes.tsv.gz` - Cell barcodes
- `features.tsv.gz` - Gene information  
- `matrix.mtx.gz` - Count matrix

### Other Formats
- `.h5` files (HDF5 format)
- `.csv` or `.tsv` files with count matrices
- `.rds` files with pre-processed R objects

## Example Data Structure

```
data/
├── experiment1/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── experiment2/
│   └── data.h5
└── metadata/
    ├── sample_info.csv
    └── experimental_design.csv
```

## Usage

Place your data files in this directory and modify the analysis scripts in the `R/` directory to point to your specific data files.

**Note:** Large data files should be stored externally (e.g., on a server or cloud storage) and downloaded as needed. Consider adding large data files to `.gitignore` to avoid committing them to version control.