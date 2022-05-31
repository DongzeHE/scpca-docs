# Single-cell gene expression file contents

Single-cell or single-nuclei gene expression data (filtered or unfiltered) is provided for use with R as an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
This object contains the expression data, cell and gene metrics, associated metadata, and, in the case of multimodal data like CITE-seq, data from additional cell-based assays.

We highly encourage you to familiarize yourself with the general object structure and functions available as part of the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor.
Below we present some details about the specific contents of the objects we provide.

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
sce <- readRDS("SCPCL000000_filtered.rds")
```

## Components of a `SingleCellExperiment` object

### Expression counts

The `counts` assay of the `SingleCellExperiment` object for single-cell and single-nuclei experiments (both unfiltered and filtered) contains the primary RNA-seq expression data as integer counts.
The data is stored as a sparse matrix, and each column represents a cell or droplet, each row a gene.
Column names are cell barcode sequences and row names are Ensembl gene IDs.
The `counts` assay can be accessed with the following R code:

```r
count_matrix <- counts(sce)
```

### Cell metrics

Cell metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `colData` slot, with the cell barcodes as the names of the rows.

```r
cell_metrics <- colData(sce)
```

The following per-cell data columns are included for each cell, calculated using the [`scuttle::addPerCellQCMetrics`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) function.

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sum`                   | UMI count for RNA-seq data                                                                                                                                                                    |
| `detected`              | Number of genes detected (gene count > 0 )                                                                                                                                                    |
| `subsets_mito_sum`      | UMI count of mitochondrial genes                                                                                                                                                              |
| `subsets_mito_detected` | Number of mitochondrial genes detected                                                                                                                                                        |
| `subsets_mito_percent`  | Percent of all UMI counts assigned to mitochondrial genes                                                                                                                                     |
| `total`                 | Total UMI count for RNA-seq data and any alternative experiments (i.e., CITE-seq)                                                                                                             |

The following are additional per-cell data columns included only in `filtered` objects.
These metrics were calculated by using [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html), a package that jointly models proportion of reads belonging to mitochondrial genes and number of unique genes detected to predict low-quality cells.

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `prob_compromised`      | Probability that a cell is compromised (i.e., dead or damaged), as calculated by `miQC`                                                                                                       |
| `miQC_pass`             | Indicates whether the cell passed the default miQC filtering. `TRUE` is assigned to cells with a low probability of being compromised (`prob_compromised` < 0.75) or [sufficiently low mitochondrial content](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html#preventing-exclusion-of-low-mito-cells).  |

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
| `detected`    | Percent of cells in which the gene was detected (gene count > 0 ) |

### Experiment metadata

Metadata associated with {ref}`data processing <processing_information:Processing information>` is included in the `metadata` slot as a list.

```r
expt_metadata <- metadata(sce)
```

| Item name           | Contents                                                                                                                       |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `salmon_version`    | Version of `salmon` used for initial mapping                                                                                   |
| `reference_index`   | Transcriptome reference file used for mapping                                                                                  |
| `total_reads`       | Total number of reads processed by `salmon`                                                                                    |
| `mapped_reads`      | Number of reads successfully mapped                                                                                            |
| `mapping_tool`      | Pipeline used for mapping and quantification (`alevin-fry` for all current data in ScPCA)                                      |
| `alevinfry_version` | Version of `alevin-fry` used for mapping and quantification                                                                    |
| `af_permit_type`    | `alevin-fry generate-permit-list` method used for filtering cell barcodes                                                      |
| `af_resolution`     | `alevin-fry quant` resolution mode used                                                                                        |
| `usa_mode`          | Boolean indicating whether quantification was done using `alevin-fry` USA mode                                                 |
| `af_num_cells`      | Number of cells reported by `alevin-fry`                                                                                       |
| `tech_version`      | A string indicating the technology and version used for the single-cell library, such as 10Xv2, 10Xv3, or 10Xv3.1              |
| `transcript_type`   | Transcripts included in gene counts: `spliced` for single-cell samples and `unspliced` for single-nuclei                       |
| `miQC_model`        | The model object that `miQC` fit to the data and was used to calculate `prob_compromised`. Only present for `filtered` objects |
| `filtering_method`  | The method used for cell filtering. One of `emptyDrops`, `emptyDropsCellRanger`, or `UMI cutoff`. Only present for `filtered` objects |
| `umi_cutoff`        | The minimum UMI count per cell used as a threshold for filtering. Only present for `filtered` objects where the `filtering_method` is `UMI cutoff` |

## Additional `SingleCellExperiment` components for CITE-seq libraries

CITE-seq data, when present, is included within the `SingleCellExperiment` as an "Alternative Experiment" named `"CITEseq"` , which can be accessed with the following command:

```r
altExp(sce, "CITEseq")
```

Within this, the main expression matrix is again found in the `counts` assay, with each column corresponding to a cell or droplet (in the same order as the parent `SingleCellExperiment`) and each row corresponding to an antibody derived tag (ADT).
Column names are again cell barcode sequences and row names the antibody targets for each ADT.

The following additional per-cell data columns for the CITE-seq data can be found in the main `colData` data frame (accessed with `colData(sce)` [as above](#cell-metrics)).

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `altexps_CITEseq_sum  `    | UMI count for CITE-seq ADTs                       |
| `altexps_CITEseq_detected` | Number of ADTs detected per cell (ADT count > 0 ) |
| `altexps_CITEseq_percent`  | Percent of `total` UMI count from ADT reads       |

Metrics for each of the ADTs assayed can be found as a `DataFrame` stored as `rowData` within the alternative experiment:

```r
adt_info <- rowData(altExp(sce, "CITEseq"))
```

This data frame contains the following columns with statistics for each ADT:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `mean`      | Mean ADT count across all cells/droplets                       |
| `detected`  | Percent of cells in which the ADT was detected (ADT count > 0 ) |

Finally, additional metadata for the CITE-seq data processing can be found in the metadata slot of the alternative experiment, with the same contents as the [parent experiment metadata](#experiment-metadata)

```r
citeseq_metadata <- metadata(altExp(sce, "CITEseq"))
```

## Additional `SingleCellExperiment` components for multiplexed libraries

Multiplexed libraries will contain a number of additional components and fields.

Hashtag oligo (HTO) quantification for each cell is included within the `SingleCellExperiment` as an "Alternative Experiment" named `"cellhash"` , which can be accessed with the following command:

```r
altExp(sce, "cellhash")
```

Within this, the main data matrix is again found in the `counts` assay, with each column corresponding to a cell or droplet (in the same order as the parent `SingleCellExperiment`) and each row corresponding to a hashtag oligo (HTO).
Column names are again cell barcode sequences and row names the HTO ids for all assayed HTOs. 

The following additional per-cell data columns for the cellhash data can be found in the main `colData` data frame (accessed with `colData(sce)` [as above](#cell-metrics)).

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `altexps_cellhash_sum  `    | UMI count for cellhash HTOs                       |
| `altexps_cellhash_detected` | Number of HTOs detected per cell (HTO count > 0 ) |
| `altexps_cellhash_percent`  | Percent of `total` UMI count from HTO reads       |

Metrics for each of the HTOs assayed can be found as a `DataFrame` stored as `rowData` within the alternative experiment:

```r
hto_info <- rowData(altExp(sce, "cellhash"))
```

This data frame contains the following columns with statistics for each HTO:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `mean`      | Mean HTO count across all cells/droplets                       |
| `detected`  | Percent of cells in which the HTO was detected (HTO count > 0 ) |
| `sample_id` | Sample ID for this library that corresponds to the HTO (only present in `_filtered.rds` files) |

Note that in the unfiltered `SingleCellExperiment` objects, this may include hashtag oligos that do not correspond to any included sample, but were part of the reference set used for mapping.

### Demultiplexing results

Demultiplexing results are included only in the `_filtered.rds` files. 
The demultiplexing methods applied for these files are as described in the {ref}`multiplex data processing section <processing_information:Multiplexed libraries>`.

Demultiplexing analysis adds the following additional fields to the `colData(sce)` data frame:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `hashedDrops_sampleid`  | Most likely sample as called be `DropletUtils::hashedDrops` 
| `HTODemux_sampleid`  | Most likely sample as called be `Seurat::HTODemux`
| `vireo_sampleid`  | Most likely sample as called by `vireo` (genetic demultiplexing)        

### Additional demultiplexing statistics

Each demultiplexing method generates additional statistics specific to the method that you may wish to access, including probabilities, alternative calls, and potential doublet information.

For methods that rely on the HTO data, these statistics are found in the `colData(altExp(sce, "cellhash"))` data frame; 
`DropletUtils::hashedDrops` statistics have the prefix `hashedDrops_` and `Seurat::HTODemux` statistics have the prefix `HTODemux`.

Genetic demultiplexing statistics are found in the main `colData(sce)` data frame, with the prefix `vireo_`.


