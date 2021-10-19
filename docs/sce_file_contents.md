# Single-cell gene expression file contents

Single cell or single-nuclei gene expression data (filtered or unfiltered) is provided for use with R as an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
This object contains the expression data, cell and gene metrics, associated metadata, and, in the case of multimodal data like CITE-seq, data from additional cell-based assays.

We highly encourage you to familiarize yourself with the general object structure and functions available as part of the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor.
Below we present some details about the specific contents of the objects we provide.

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
sce <- readRDS("SCPCL000000_filtered.rds")
```

## Common fields

### Expression counts

For both unfiltered and filtered single-cell and single-nucleotide experiments, the primary cell by gene RNA-seq expression count data is provided as the `"counts"` assay.
Each column represents a cell or droplet, each row a gene.
Column names are cell barcode sequences and row names are Ensembl gene IDs. 
This data is stored as a sparse matrix that can be accessed with the following R code:

```r
count_matrix <- counts(sce)
```

### Cell metrics

Cell metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `colData` slot, with the cell barcodes as the names of the rows.

```r
cell_metrics <- colData(sce)
```

The following data columns are included for all cells, calculated using the [`scuttle::addPerCellQCMetrics`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) function. 

| Column name             | Contents                                                                          |
| ----------------------- | --------------------------------------------------------------------------------- |
| `sum`                   | Total UMI count for RNA-seq data                                                  |
| `detected`              | Number of genes detected per cell (gene count > 0 )                               |
| `subsets_mito_sum`      | UMI count of mitochondrial genes                                                  |
| `subsets_mito_detected` | Number of mitochondrial genes detected                                            |
| `subsets_mito_percent`  | Percent of all UMI counts assigned to mitochondrial genes                         |
| `total`                 | Total UMI count for RNA-seq data and any alternative experiments (i.e., CITE-seq) |


### Gene information and metrics

Gene information and metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `rowData` slot, with the Ensembl ID as the names of the rows.

```r
gene_info <- rowData(sce)
```

The following columns are included for all genes. 
Metrics were calculated using the [`scuttle::addPerFeatureQCMetrics`](https://rdrr.io/github/LTLA/scuttle/man/addPerFeatureQCMetrics.html) function.

| Column name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `gene_symbol` | [HUGO](https://www.genenames.org) gene symbol, if defined        |
| `mean`        | Mean count across all cells/droplets                             |
| `detected`    | Number of cells in which the gene was detected (gene count > 0 ) |

### Experiment metadata

Metadata associated with data processing is included in the `metadata` slot as a list
