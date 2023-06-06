# Getting started with an ScPCA dataset

This section provides information on next steps you might take after downloading a dataset from the ScPCA portal.

We also have a separate GitHub repository that contains workflows for some common analysis performed on single-cell RNA-sequencing data.
These workflows are designed to apply the same analysis (e.g., clustering) across multiple samples in parallel.
These workflows and more resources for processing single-cell and single-nuclei datasets can be found in the [`scpca-downstream-analyses` repository](https://github.com/AlexsLemonade/scpca-downstream-analyses).

## Importing ScPCA data into R

Quantified single-cell or single-nuclei gene expression data is provided as an RDS file as described in the {ref}`single cell gene expression file contents section<sce_file_contents:single-cell gene expression file contents>`.
There are three RDS files that are available for each library: an `unfiltered.rds`, a `filtered.rds`, and a `processed.rds` file.
The `unfiltered.rds` file contains the gene expression data for all droplets, regardless of the presence of a cell or not.
The `filtered.rds` files contains the gene expression data for only droplets that are likely to contain cells, removing any probable empty droplets.
See the {ref}`section on filtering cells<processing_information:filtering cells>` for more information on how we remove potential empty droplets.

The `processed.rds` files are further filtered to remove any low quality cells and contain both the raw and normalized gene expression data for the identified cells.
See the description of the {ref}`processed gene expression data <processing_information:Processed gene expression data>` for more information on the `processed` objects.

In most scenarios, we recommend starting with the `processed.rds` file.

The first step in analyzing the provided gene expression data will be to import the data into R.
As a reminder, each RDS file contains a [`SingleCellExperiment` object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).
Refer to {ref}`single-cell gene expression file contents<sce_file_contents:single-cell gene expression file contents>` for a more detailed description on the contents of the included `SingleCellExperiment` object.

To work with `SingleCellExperiment` objects in R, we need to ensure that we have the [`SingleCellExperiment` package](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) installed and loaded.

The following commands can be used to import the RDS file into R and save the `SingleCellExperiment` object:

```r
# if SingleCellExperiment is not installed, install the package
# otherwise the installation step can be skipped
if (!("SingleCellExperiment" %in% installed.packages())) {
  BiocManager::install("SingleCellExperiment")
}

library(SingleCellExperiment)
# read in the RDS file, including the path to the file's location
processed_sce <- readRDS("SCPCS000000/SCPCL000000_processed.rds")
```

More resources for learning about `SingleCellExperiment` objects:

- [`SingleCellExperiment` objects Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) on working with `SingleCellExperiment` objects
- [Orchestrating Single Cell Analysis chapter on the `SingleCellExperiment` class](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html)

## Working with the processed `SingleCellExperiment` objects

The `SingleCellExperiment` objects stored in the `processed.rds` files have undergone additional quality control steps to remove low quality cells.
In addition, normalized expression counts and dimensionality reduction (principal component analysis and UMAP) have been calculated.

The following commands can be used to access the raw and normalized count matrices:

```r
# the raw counts matrix stored in the processed object
raw_counts <- counts(processed_sce)

# log normalized counts matrix stored in the processed object
normalized_counts <- logcounts(processed_sce)
```

If you have ADT data from a CITE-seq experiment, you can use the following commands to access raw and normalized ADT count matrices:

```r
# the raw ADT counts matrix stored in the processed object
raw_adt_counts <- counts(altExp(processed_sce))

# log normalized ADT counts matrix stored in the processed object
# if this value is NULL, it means normalization failed
normalized_adt_counts <- logcounts(altExp(processed_sce))
```


Dimensionality reduction results can be accessed using the following command:

```r
# extract principal component results
pca_results <- reducedDim(processed_sce, "PCA")

# extract UMAP results
umap_results <- reducedDim(processed_sce, "UMAP")
```

Principal components were calculated from a set of highly variable genes identified for a given library.
The list of highly variable genes used for this calculation, in order from highest to lowest variation, is stored in the `metadata` of the `SingleCellExperiment` object.

This list can be accessed using the following command:

```r
highly_variable_genes <- metadata(processed_sce)$highly_variable_genes
```

This data is immediately ready for clustering and further analysis to answer biological questions of interest.

See these resources for more information on clustering:
 - [Clustering chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.14/OSCA.basic/clustering.html)
 - [Quantifying clustering behavior in Orchestrating Single Cell Analysis](https://bioconductor.org/books/release/OSCA.advanced/clustering-redux.html#quantifying-clustering-behavior)

## Working with the filtered `SingleCellExperiment` objects

If you prefer to work with the `filtered` objects and perform the quality control and normalization yourself, see below for details on working with these objects.

### Quality control

Before performing any downstream steps, we recommend removing low quality cells from the dataset.
Low quality cells include those with a higher percentage of reads from mitochondrial genes (i.e., those that are damaged or dying) and those with a lower number of total reads and unique genes identified (i.e., those with inefficient reverse transcription or PCR amplification).

All `filtered.rds` objects include [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html) results found in the {ref}`colData() of the SingleCellExperiment object<sce_file_contents:cell metrics>`.
`miQC` jointly models the proportion of mitochondrial reads and the number of unique genes detected in each cell to calculate the probability of a cell being compromised (i.e., dead or damaged).
High-quality cells are those with a low probability of being being compromised (< 0.75) or sufficiently low mitochondrial content.
All cells that are identified as low-quality cells will have `FALSE` in the `miQC_pass` column of the `colData()` and can be removed prior to downstream analyses.

The following command can be used to remove the low-quality cells:

```r
# read in the filtered object
sce <- readRDS("SCPCS000000/SCPCL000000_filtered.rds")

# filter the `SingleCellExperiment`
# the `$` notation denotes access to the colData slot of the `SingleCellExperiment` object
filtered_sce <- sce[, which(sce$miQC_pass)]
```

For more information on using `miQC` for filtering cells, see the following resources:
- [`miQC` vignette](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html)
- [Hippen _et al._ 2021](https://doi.org/10.1371/journal.pcbi.1009290)

You can also directly filter cells based on the number of unique genes, total reads, and fraction of mitochondrial reads.
We used the function [`scuttle::addPerCellQCMetrics()`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) to calculate these metrics and have included them in the `colData` of the filtered `SingleCellExperiment` objects.

The following columns are added by `scuttle::addPerCellQCMetrics()` and can be found in the `colData`:

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------|
| `sum`                   | UMI count for RNA-seq data                                                                                                                           |
| `detected`              | Number of genes detected (gene count > 0 )                                                                                                                                                    |
| `subsets_mito_sum`      | UMI count of mitochondrial genes                                                                                                                                                              |
| `subsets_mito_detected` | Number of mitochondrial genes detected                                                                                                                                                        |
| `subsets_mito_percent`  | Percent of all UMI counts assigned to mitochondrial genes                                                                                                                                     |
| `total`                 | Total UMI count for RNA-seq data and any alternative experiments (i.e., ADT data from CITE-seq)                                                                                                             |

These metrics can be used to directly filter the `SingleCellExperiment` object based on informed thresholds.
If you are planning to filter low quality cells using such thresholds, we encourage you to read more about the various metrics and plot the distribution of each metric before deciding on which cells to exclude.
The [Quality Control chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/quality-control.html#quality-control) provides a nice guide to checking diagnostic plots and then choosing cutoffs.

If you have ADT data from a CITE-seq experiment, please refer to the section [Special considerations for CITE-seq experiments](#special-considerations-for-cite-seq-experiments) below for information on how to filter cells based on ADT-level statistics.
### Normalization

The provided data contains unnormalized raw counts.
We recommend using the `scran` and `scater` packages to add normalized counts to the `SingleCellExperiment` object.

Follow these steps to perform normalization using `scran` and `scater`:

```r
# Cluster similar cells
qclust <- scran::quickCluster(filtered_sce)

# Compute sum factors for each cell cluster grouping.
filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)

# Normalize and log transform.
normalized_sce <- scater::logNormCounts(filtered_sce)
```

Here we provide more resources on understanding normalization in single-cell RNA-seq analysis:

- [Normalization chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.14/OSCA.basic/normalization.html)
- [Hemberg lab scRNA-seq course section on normalization methods](https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html#normalization-theory)
- [Stegle _et al._ (2015) Computational and analytical challenges in single-cell transcriptomics](https://doi.org/10.1038/nrg3833).  Includes a discussion of normalization and technical variance in scRNA-seq.
- [Lun _et al._ (2016) Pooling across cells to normalize single-cell RNA sequencing data with many zero counts](https://doi.org/10.1186/s13059-016-0947-7)

If you have ADT data from a CITE-seq experiment, please refer to the section [Special considerations for CITE-seq experiments](#special-considerations-for-cite-seq-experiments) below for information on how to normalize ADT counts.

### Dimensionality Reduction

Dimensionality reduction is commonly used as a precursor to plotting, clustering, and other downstream analysis.

It is common to start with performing [principal component analysis (PCA)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#principal-components-analysis), a technique that identifies new axes that capture the largest amount of variation in the data.
The PCA results can be calculated and stored in the `SingleCellExperiment` object using the following command:

```r
# calculate a PCA matrix using the top 500 most highly variable genes (default)
normalized_sce <- scater::runPCA(normalized_sce, ntop = 500)
```

Here we are calculating PCA by using the default of the top 500 most highly variable genes as input, however this is not always the optimal choice.
We encourage you to visit the [Feature selection chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/feature-selection.html#feature-selection) to read more about modeling gene variance and selecting the highly variable genes.

PCA is commonly used for initial dimensionality reduction.
However, we can use more advanced techniques, like [UMAP (Uniform Manifold Approximation and Projection)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#uniform-manifold-approximation-and-projection), that may be better for visualization.
Use caution when interpreting UMAP results, as location, spacing, and density of clusters can be dependent on parameter choices and random effects and do not always accurately represent the relationships among cells.

UMAP can also be quite slow for a large dataset, so we can use the previous PCA results as input to speed up the analysis.

```r
# Run UMAP using already stored PCA results
normalized_sce <- scater::runUMAP(normalized_sce, dimred = "PCA")
```

See below for more resources on dimensionality reduction:

- [Dimensionality Reduction chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#dimensionality-reduction)
- [Dimension Reductions in R](https://rpubs.com/Saskia/520216)
- [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

## What if I want to use Seurat?

The files available for download contain [`SingleCellExperiment` objects](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
If desired, these can be converted into Seurat objects.
You can find the code needed to convert the `SingleCellExperiment` objects to Seurat objects in the {ref}`FAQ section<FAQ:what if i want to use seurat instead of bioconductor?>`.

After converting the object to a Seurat object, the same steps outlined above (quality control, filtering, normalization, dimensionality reduction) can be followed but using functions available as part of the Seurat package.

Here are some resources that can be used to get you started working with Seurat objects:
- [Getting started tutorial in Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) including quality control, normalization, and dimensionality reduction.
- [Converting Seurat objects to and from a `SingleCellExperiment`](https://satijalab.org/seurat/articles/conversion_vignette.html)
- [Harvard Chan Bioinformatics Core course on single-cell RNA-seq analysis with Seurat](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html)


## What if I want to use scanpy?

There are a variety of ways to convert the counts data into a python-compatible format.
We have found that one of the more efficient ways is to first convert the counts data into 10x format using [`DropletUtils::write10xCounts()`](https://rdrr.io/bioc/DropletUtils/man/write10xCounts.html) and then read the files into Python using the [`scanpy` package](https://scanpy.readthedocs.io/en/stable/).

You can find the code needed to perform these steps outlined in the {ref}`FAQ section on using Python<FAQ:what if i want to use python instead of r?>`.

After creating an `AnnData` object in `scanpy`, the same steps outlined above (quality control, filtering, normalization, dimensionality reduction) can be followed but using functions available as part of the `scanpy` package.

Here are some resources that can be used to get you started working with `AnnData` objects:
- [Preprocessing and clustering tutorial in scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
- [Converting directly from R to Python](https://theislab.github.io/scanpy-in-R/#converting-from-r-to-python)
- [Homepage for scanpy tutorials](https://scanpy.readthedocs.io/en/latest/tutorials.html)


## Special considerations for CITE-seq experiments

### Filtering cells based on ADT quality control

The `SingleCellExperiment` objects stored in both the `filtered.rds` and the `processed.rds` have not been filtered based on ADT-level statistics, but both objects contain indicator variables to quickly perform this filtering.

In both objects, you can see cells flagged for removal based on ADT-level QC statistics in the alternative experiment's `discard` column, where `TRUE` indicates the cell should be removed.

```r
# See which cells should be removed (`TRUE`) vs. kept (`FALSE`)
# This information is in both `filtered.rds` and `processed.rds`
altExp(sce)$discard
```

In the `processed.rds` object specifically, this information is also recorded as values "Remove" and "Keep":

```r
# See which cells should be removed vs. kept in `processed.rds`
processed_sce$adt_scpca_filter
```
### Normalizing ADT counts

The `filtered.rds` object contains a raw ADT counts matrix, but not a matrix of log normalized ADT expression.
To perform normalization yourself, we first recommend filtering cells flagged for removal by ADT-level QC statistics to reduce the chances of normalization failing:

```r
# Filter cells based on ADT-level QC statistics to
# prepare for normalization
filtered_sce <- filtered_sce[, adt_scpca_filter == "Keep"]
```

Normalization can then be performed with `scuttle` and `scater`:

```r
# Calculate size factors for use in ADT normalization,
# using the previously-calculated ambient profile
altExp(filtered_sce) <- scuttle::computeMedianFactors(
  altExp(filtered_sce),
  reference = metadata(altExp(filtered_sce))$ambient_profile
)

# Normalize and log transform
altExp(filtered_sce) <- scater::logNormCounts(altExp(filtered_sce))
```

A normalized expression matrix is also provided in the `processed.rds` object, specifically for all cells which are indicated as `"Keep"` in the `processed_sce$adt_scpca_filter` column:

```r
# Access log normalized ADT expression matrix
logcounts(altExp(processed_sce))
```

Here are some additional resources that can be used for working with ADT counts from CITE-seq experiments:
- [Integrating with Protein Abundance, Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html)
- [Seurat vignette on using with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette.html)



## Special considerations for multiplexed samples

If the dataset that you have downloaded contains samples that were multiplexed (i.e. cells from multiple samples have been combined together into one library), then the steps to take after downloading will be a little different than for libraries that only contain cells corresponding to a single sample.
Here, multiplexed samples refer to samples that have been combined together into one library using cell hashing and then sequenced together.
This means that a single library contains cells or nuclei that correspond to multiple samples.
Each sample has been tagged with a hashtag oligo (HTO) prior to mixing, and that HTO can be used to identify which cells belong to which sample within a multiplexed library.
The libraries available for download on the portal have not been separated by sample (i.e. demultiplexed), and therefore contain data from multiple samples.

Libraries containing multiplexed samples can be initially processed using the same workflow described above including removal of [low quality cells](#quality-control), [normalization](#normalization), and [dimensionality reduction](#dimensionality-reduction).
Demultiplexing can then be used to identify the sample that each cell is from.
Demultiplexing has already been performed using both [`Seurat::HTODemux`](https://satijalab.org/seurat/reference/htodemux) and [`DropletUtils::hashedDrops`](https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html).
For samples where corresponding bulk RNA-sequencing data is available, {ref}`genetic demultiplexing <processing_information:Genetic demultiplexing>` was also conducted.
The results from demultiplexing using these methods have been summarized and are present in the `colData` of the `SingleCellExperiment` object in the `_filtered.rds` file only.
The `hashedDrops_sampleid`, `HTODemux_sampleid`, and `vireo_sampleid` columns in the `colData` report the sample called for each cell by the specified demultiplexing method.
If a confident call was not made for a cell by the demultiplexing method, the column will have a value of `NA`.
For more information on how to access the full demultiplexing results, see {ref}`this description of demultiplexing results <sce_file_contents:demultiplexing results>`.

If desired, the sample calls identified from demultiplexing can be used to separate the `SingleCellExperiment` object by sample for downstream analysis.

```r
# view the summary of the sample calls for genetic demultiplexing
# genetic demultiplexing results are stored in the vireo_sampleid column found in the colData
table(multiplexed_sce$vireo_sampleid)

# identify cells that are assigned to sample A (use `which` to remove NA)
sampleA_cells <- which(multiplexed_sce$vireo_sampleid == "sampleA")

# create a new sce that only contains cells from sample A
sampleA_sce <- multiplexed_sce[, sampleA_cells]
```

Here are some additional resources that can be used for working with multiplexed samples (or those with cell hashing):
- [Demultiplexing on HTO Abundance, Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.14/OSCA.advanced/droplet-processing.html#demultiplexing-on-hto-abundance)
- [Seurat vignette on demultiplexing with Seurat::HTODemux](https://satijalab.org/seurat/articles/hashing_vignette.html)
- [Weber _et al._ (2021) Genetic demultiplexing of pooled single-cell RNA-sequencing samples in cancer facilitates effective experimental design](https://doi.org/10.1093/gigascience/giab062)
