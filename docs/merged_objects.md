# Merged objects

Each merged object contain _combined information_ for all individual libraries in a given ScPCA project.
While each individual object, as described on the {ref}`Single-cell gene expression file contents page <sce_file_contents>`, contains quantified gene expression results for a single library, each merged object contains all gene expression results, including gene expression counts and metadata, for all libraries in a given ScPCA project
This information includes quantified gene expression data, cell and gene metrics, and associated metadata for all libraries.
See {ref}`the section on merged object processing <processing_information:merged objects` for more information on how these objects were prepared.


Merged objects are provided in two formats:
  - As an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.17/OSCA.intro/the-singlecellexperiment-class.html) for use in R.
  - As an HDF5 file containing an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html) for use in Python.


Below we present some details about the specific contents of the objects we provide.

## Components of a SingleCellExperiment merged object

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
merged_sce <- readRDS("SCPCP000000_merged.rds")
```

### SingleCellExperiment expression counts

Merged `SingleCellExperiment` objects contain two main assays, `counts` and `logcounts`, each containing RNA-seq expression data for all libraries in a given ScPCA project combined into a single matrix.
The `counts` assay contains the primary raw counts represented as integers, and the `logcounts` assay contains normalized counts as described in {ref}`the data post-processing section <processing_information:processed gene expression data>`.

The `counts` assay includes reads aligned to both spliced and unspliced cDNA (see the section on {ref}`Post Alevin-fry processing <processing_information:post alevin-fry processing>`).
Each assay is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `counts` and `logcounts` assays can be accessed with the following R code:

```r
counts(merged_sce) # combined counts matrix
logcounts(merged_sce) # combined logcounts matrix
```

Column names are cell barcode sequences prefixed with the originating library id, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed with the following R code:

```r
colnames(merged_sce) # matrix column names
rownames(merged_sce) # matrix row names
```

There is also a `spliced` assay which contains the counts matrix with only reads from spliced cDNA.

### SingleCellExperiment cell metrics


Cell metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `colData` slot, where row names are the cell barcode prefixed with the originating library id, e.g. `SCPCL000000-{barcode}`.

```r
colData(merged_sce) # cell metrics
```

The following per-cell data columns are included for each cell, calculated using the [`scuttle::addPerCellQCMetrics()`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) function.

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sum`                   | UMI count for RNA-seq data                                                                                                                                                                    |
| `detected`              | Number of genes detected (gene count > 0 )                                                                                                                                                    |
| `subsets_mito_sum`      | UMI count of mitochondrial genes                                                                                                                                                              |
| `subsets_mito_detected` | Number of mitochondrial genes detected                                                                                                                                                        |
| `subsets_mito_percent`  | Percent of all UMI counts assigned to mitochondrial genes                                                                                                                                     |
| `total`                 | Total UMI count for RNA-seq data and any alternative experiments (i.e., ADT data from CITE-seq)                                                                                                             |
| `prob_compromised`      | Probability that a cell is compromised (i.e., dead or damaged), as calculated by `miQC`                                                                                                       |
| `miQC_pass`             | Indicates whether the cell passed the default miQC filtering. `TRUE` is assigned to cells with a low probability of being compromised (`prob_compromised` < 0.75) or [sufficiently low mitochondrial content](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html#preventing-exclusion-of-low-mito-cells)  |
| `scpca_filter` | Labels cells as either `Keep` or `Remove` based on filtering criteria (`prob_compromised` < 0.75 and number of unique genes detected > 200) |
| `adt_scpca_filter` | If CITE-seq was performed, labels cells as either `Keep` or `Remove` based on ADT filtering criteria (`discard = TRUE` as determined by [`DropletUtils::CleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html)) |
| `submitter_celltype_annotation` | If available, cell type annotations obtained from the group that submitted the original data. Cells that the submitter did not annotate are labeled as `"Submitter-excluded"` |
| `singler_celltype_annotation`  | If cell typing with `SingleR` was performed, the annotated cell type. Cells labeled as `NA` are those which `SingleR` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"` |
| `singler_celltype_ontology`  | If cell typing with `SingleR` was performed with ontology labels, the annotated cell type's ontology ID. Cells labeled as `NA` are those which `SingleR` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"` |
| `cellassign_celltype_annotation`  | If cell typing with `CellAssign` was performed, the annotated cell type. Cells labeled as `"other"` are those which `CellAssign` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"`  |
| `cellassign_max_prediction`  | If cell typing with `CellAssign` was performed, the annotation's prediction score (probability)  |





## Components of an AnnData merged object

Before getting started, we highly encourage you to familiarize yourself with the general `AnnData` object structure and functions available as part of the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html).
For the most part, the `AnnData` objects that we provide are formatted to match the expected data format for [`CELLxGENE`](https://cellxgene.cziscience.com/) following [schema version `3.0.0`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md).

To begin, you will need to load the `AnnData` package and read the HDF5 file:

```python
import anndata
merged_adata_object = anndata.read_h5ad("SCPCP000000_merged_rna.hdf5")
```

### AnnData expression counts

Merged `AnnData` objects contain two data matrices, each containing RNA-seq expression data for all libraries in a given ScPCA project combined into a single matrix.
The data matrix `raw.X` of the merged `AnnData` object contains the RNA-seq expression data as primary integer counts, and the data matrix `X` contains the RNA-seq expression data as normalized counts.
The data is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `raw.X` and `X` matrices can be accessed with the following python code:

```python
merged_adata_object.raw.X # raw count matrix
merged_adata_object.X # normalized count matrix
```

Column names are cell barcode sequences prefixed with the originating library id, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed as with the following python code:


```python
merged_adata_object.obs_names # matrix column names
merged_adata_object.var_names # matrix row names
```


