# Merged objects

These objects contain processed gene expression data, cell and gene metrics, associated metadata, and, in the case of multimodal data like ADTs from CITE-seq experiments, data from additional cell-based assays for all libraries in a given ScPCA project.
See {ref}`the section on merged object processing <processing_information:merged objects` for more information on how these objects were prepared.

Merged objects are provided in two formats:
  - As an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.17/OSCA.intro/the-singlecellexperiment-class.html) for use in R.
  - As an HDF5 file containing an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html) for use in Python.


Below we present some details about the specific contents of the objects we provide.

## Components of a SingleCellExperiment merged object

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
merged_sce <- readRDS("SCPCL000000_merged.rds")
```

### SingleCellExperiment expression counts

The `counts` and `logcounts` assays of the `SingleCellExperiment` object for single-cell and single-nuclei experiments contains the primary RNA-seq expression data.
The `counts` assay contains raw counts represented as integers, and the `logcounts` assay contains normalized (on a per-library basis) counts as described in {ref}`the data post-processing section <processing_information:processed gene expression data>`.

Both assays include reads aligned to both spliced and unspliced cDNA (see the section on {ref}`Post Alevin-fry processing <processing_information:post alevin-fry processing>`).
The data is stored as a sparse matrix, and each column represents a cell or droplet, and each row represents a gene.
Column names are cell barcode sequences prefixed with the originating library id, e.g. `SCPCL00000-{barcode}`, and row names are Ensembl gene IDs.
The `counts` and `logcounts` assays can be accessed with the following R code:

```r
counts(merged_sce) # counts matrix
logcounts(merged_sce) # logcounts matrix
```

Additionally, the `spliced` assay contains a counts matrix that includes reads from spliced cDNA only.


