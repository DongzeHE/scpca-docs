# Getting started with an ScPCA dataset

This section provides information on next steps you might take after downloading a dataset from the ScPCA portal.
Quantified single-cell or single-nuclei gene expression data is provided as either [`SingleCellExperiment` objects (`.rds` files)](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) or [`AnnData` objects (`.hdf5` files)](https://anndata.readthedocs.io/en/latest/index.html).
A full description of the contents of the `SingleCellExperiment` and `AnnData` objects can be found in the {ref}`single cell gene expression file contents section<sce_file_contents:single-cell gene expression file contents>`.

There are three objects available for each library: an unfiltered object (`unfiltered.rds` or `unfiltered_rna.hdf5`), a filtered object (`filtered.rds` or `filtered_rna.hdf5`), and a processed object (`processed.rds` or `processed_rna.hdf5`).
The unfiltered object contains the gene expression data for all droplets, regardless of the presence of a cell or not.
The filtered object contains the gene expression data for only droplets that are likely to contain cells, removing any probable empty droplets.
See the {ref}`section on filtering cells<processing_information:Filtering cells>` for more information on how we remove potential empty droplets.

The processed objects are further filtered to remove any low quality cells and contain both the raw and normalized gene expression data for the identified cells.
Additionally, these objects contain dimensionality reduction results (principal component analysis and UMAP) and clustering assignments.
See the description of the {ref}`processed gene expression data <processing_information:Processed gene expression data>` for more information on the `processed` objects.

We recommend starting with the processed objects.
  - For working with the processed `SingleCellExperiment`object found in the `processed.rds` file see [Importing ScPCA data from `SingleCellExperiment` objects into R below](#importing-scpca-data-from-singlecellexperiment-objects-into-r).
  - For working with the processed `AnnData` objects found in the `processed_rna.hdf5` file see [Importing ScPCA data from `AnnData` objects into Python below](#importing-scpca-data-from-anndata-objects-into-python).

Note: We also have a separate GitHub repository that contains workflows for some common analysis performed on single-cell RNA-sequencing data.
These workflows are designed to apply the same analysis (e.g., clustering) across multiple samples in parallel and is currently only for use with `SingleCellExperiment` objects.
These workflows and more resources for processing single-cell and single-nuclei datasets can be found in the [`scpca-downstream-analyses` repository](https://github.com/AlexsLemonade/scpca-downstream-analyses).

## Importing ScPCA data from `SingleCellExperiment` objects into R

The first step in analyzing the provided gene expression data stored in the `SingleCellExperiment` objects will be to import the data into R.
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

## Importing ScPCA data from `AnnData` objects into Python

The first step in analyzing the provided gene expression data stored in the `AnnData` objects will be to import the data into python.
To work with `AnnData` objects in python, we need to ensure that we have the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html) installed and loaded.

The following commands can be used to import the HDF5 file into python and save the `AnnData` object:

```python
# load anndata module
import anndata

# read in the HDF5 file, including the path to the file's location
processed_adata = anndata.read_h5ad("SCPCS000000/SCPCL000000_processed_rna.hdf5")
```

More resources for learning about `AnnData` objects:

- [Getting started with `AnnData` objects](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [Preprocessing and clustering tutorial in scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
- [Homepage for scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html)

## Working with processed objects

### Gene expression data

The processed objects contain both the raw and normalized gene expression data.
The following commands can be used to access the expression data found in the processed `SingleCellExperiment` objects:

```r
# the raw counts matrix stored in the processed object
counts(processed_sce)

# log normalized counts matrix stored in the processed object
logcounts(processed_sce)
```

The following commands can be used to access the expression data found in the processed `AnnData` objects:

```python
# the raw counts matrix stored in the processed object
processed_adata.raw.X

# log normalized counts matrix stored in the processed object
processed_adata.X
```

Here we provide more resources on understanding normalization in single-cell RNA-seq analysis:

- [Normalization chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.14/OSCA.basic/normalization.html)
- [Hemberg lab scRNA-seq course section on normalization methods](https://www.singlecellcourse.org/basic-quality-control-qc-and-exploration-of-scrna-seq-datasets.html#normalization-theory)
- [Stegle _et al._ (2015) Computational and analytical challenges in single-cell transcriptomics](https://doi.org/10.1038/nrg3833).  Includes a discussion of normalization and technical variance in scRNA-seq.
- [Lun _et al._ (2016) Pooling across cells to normalize single-cell RNA sequencing data with many zero counts](https://doi.org/10.1186/s13059-016-0947-7)

### Quality control data

The processed `SingleCellExperiment` objects already have undergone additional quality control steps to remove low quality cells.
Low quality cells include those with a higher percentage of reads from mitochondrial genes (i.e., those that are damaged or dying) and those with a lower number of total reads and unique genes identified (i.e., those with inefficient reverse transcription or PCR amplification).

All processed objects include [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html) results found in the {ref}`colData() of the SingleCellExperiment object<sce_file_contents:Singlecellexperiment cell metrics>` or the {ref}`.obs slot of the AnnData object<sce_file_contents:Anndata cell metrics`.
`miQC` jointly models the proportion of mitochondrial reads and the number of unique genes detected in each cell to calculate the probability of a cell being compromised (i.e., dead or damaged).
High-quality cells are those with a low probability of being being compromised (< 0.75) or sufficiently low mitochondrial content.

The probability compromised for each cell as calculated by `miQC` can be found in the `prob_compromised` column of the `colData` of `SingleCellExperiment` objects.

```r
# probability compromised for each cell
processed_sce$prob_compromised
```

The probability compromised for each cell as calculated by `miQC` can be found in the `prob_compromised` column of the `.obs` slot of `AnnData` objects.

```python
# probability compromised for each cell
processed_adata.obs["prob_compromised"]
```

Additionally we include a column, `scpca_filter`, that labels cells as either `Keep` or `Remove` based on having both a `prob_compromised` < 0.75 and number of unique genes identified > 200.
All cells included in the processed object will have `Keep` in the `scpca_filter` column.
If you prefer to work with the object prior to removal of any low-quality cells, please use the filtered object, which contains all cells that were not discarded as empty droplets.

### Dimensionality reduction

Dimensionality reduction is commonly used as a precursor to plotting, clustering, and other downstream analysis.

The processed objects contain results from performing [principal component analysis (PCA)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#principal-components-analysis), a technique that identifies new axes that capture the largest amount of variation in the data, and [uniform manifold approximation and projection(UMAP)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#uniform-manifold-approximation-and-projection), which may be better for visualization.
Use caution when interpreting UMAP results, as location, spacing, and density of clusters can be dependent on parameter choices and random effects and do not always accurately represent the relationships among cells.

Dimensionality reduction results can be accessed in the `SingleCellExperiment` objects using the following commands:

```r
# principal component analysis results
reducedDim(processed_sce, "PCA")

# UMAP results
reducedDim(processed_sce, "UMAP")
```

Dimensionality reduction results can be accessed in the `AnnData` objects using the following command:

```python
# principal component analysis results
processed_adata.obsm["X_PCA"]

# UMAP results
processed_adata.obsm["X_UMAP"]
```

See below for more resources on dimensionality reduction:

- [Dimensionality Reduction chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#dimensionality-reduction)
- [Dimension Reductions in R](https://rpubs.com/Saskia/520216)
- [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

### Highly variable genes

In the processed objects, principal components were calculated from a set of highly variable genes identified for a given library.
We encourage you to visit the [Feature selection chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/feature-selection.html#feature-selection) to read more about modeling gene variance and selecting the highly variable genes.

The list of highly variable genes used for this calculation, in order from highest to lowest variation, is stored in the `metadata` of the `SingleCellExperiment` object and the `uns` slot of the `AnnData` object.

This list can be accessed using the following command in the `SingleCellExperiment` objects:

```r
# list of highly variable genes
metadata(processed_sce)$highly_variable_genes
```

This list can be accessed using the following command in the `AnnData` objects:

```python
# list of highly variable genes
processed_adata.uns["highly_variable_genes"]
```

### Clustering

Cluster assignments obtained from [Graph-based clustering](http://bioconductor.org/books/3.16/OSCA.basic/clustering.html#clustering-graph) is also available in the processed objects.
Here clustering was performed using the Louvain algorithm with 20 nearest neighbors and Jaccard weighting.

To access the cluster assignments in the `SingleCellExperiment` object, use the following command:

```r
# cluster assignment for each cell
processed_sce$cluster
```

To access the cluster assignments in the `AnnData` object, use the following command:

```python
# cluster assignment for each cell
processed_adata.obs["cluster"]
```

See these resources for more information on clustering:
 - [Clustering chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.14/OSCA.basic/clustering.html)
 - [Quantifying clustering behavior in Orchestrating Single Cell Analysis](https://bioconductor.org/books/release/OSCA.advanced/clustering-redux.html#quantifying-clustering-behavior)


## What if I want to use Seurat?

The files available for download that contain [`SingleCellExperiment` objects](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) can also be converted into Seurat objects.
You can find the code needed to convert the `SingleCellExperiment` objects to Seurat objects in the {ref}`FAQ section<FAQ:what if i want to use seurat instead of bioconductor?>`.

After converting the object to a Seurat object, the same steps outlined above (quality control, filtering, normalization, dimensionality reduction) can be followed but using functions available as part of the Seurat package.

Here are some resources that can be used to get you started working with Seurat objects:
- [Getting started tutorial in Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) including quality control, normalization, and dimensionality reduction.
- [Converting Seurat objects to and from a `SingleCellExperiment`](https://satijalab.org/seurat/articles/conversion_vignette.html)
- [Harvard Chan Bioinformatics Core course on single-cell RNA-seq analysis with Seurat](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html)


## Special considerations for CITE-seq experiments

If the dataset you downloaded contains samples with ADT data from a CITE-seq experiment, the raw and normalized ADT expression matrices will be provided.

For `SingleCellExperiment` objects with ADT data, the ADT expression matrices will be stored as an `altExp` named `"adt"` in the same object containing the RNA expression data.

For `AnnData` objects with ADT data, the ADT expression matrices will be provided in separate files corresponding to the same three stages of data processing: an unfiltered object (`_unfiltered_adt.hdf5`), a filtered object (`_filtered_adt.hdf5`), and a processed object (`_processed_adt.hdf5`).
These files will only contain ADT expression data and not RNA expression data.

We recommend working with the processed objects as they contain both the raw and normalized ADT expression matrices along with additional quality control metrics for the CITE-seq experiment.
To access the ADT expression matrices in `SingleCellExperiment` objects, use the following commands:

```r
# View the ADT alternative experiment
# Note that the altExp name is optional, as this is the only altExp present
altExp(processed_sce, "adt")

# the raw ADT counts matrix
counts(altExp(processed_sce))

# the log normalized ADT counts matrix
logcounts(altExp(processed_sce))
```

To access the ADT matrices in `AnnData` objects you will first need to read in the `_adt.hdf5` file, and then you can access the `raw.X` and `X` matrices as shown below:

```python
import anndata
# read in the HDF5 file with ADT data, including the path to the file's location
adt_adata = anndata.read_h5ad("SCPCS000000/SCPCL000000_processed_adt.hdf5")

# the raw ADT counts matrix
adt_adata.raw.X

# the log normalized ADT counts matrix
adt_adata.X
```

Be aware that the processed objects have been filtered to remove low-quality cells based on RNA expression but have not been filtered based on ADT counts.

### Filtering cells based on ADT quality control

The `adt_scpca_filter` column indicates which cells should be removed before proceeding with downstream analyses of the ADT data, as determined by [`DropletUtils::CleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html).
This process identified cells with high levels of ambient contamination and/or high levels of negative control ADTs (if available).
Cells are labeled either as `"Keep"` (cells to retain) or `"Remove"` (cells to filter out).

To filter cells based on this column in the `SingleCellExperiment` objects, use the following command:

```r
# Filter cells based on ADT QC statistics
processed_sce[, which(processed_sce$adt_scpca_filter == "Keep")]
```

To filter cells based on this column in the `AnnData` objects, use the following command:

```python
# Filter cells based on ADT QC statistics
adt_adata[adt_adata.obs["adt_scpca_filter"] == "Keep"]
```

Note that the normalized ADT expression matrix only contains values for cells labeled as `"Keep"` in the `adt_scpca_filter` column.
Any cells labeled `"Remove"` have `NA` values in the normalized expression matrix (see {ref}`Processed ADT Data <processing_information:Processed ADT data>` for more details).

Alternatively, you can also filter cells out based on your own criteria.
For `SingleCellExperiment` objects, quality-control statistics calculated by [`DropletUtils::CleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html) are provided in the `colData` slot of the `altExp` (`colData(altExp(filtered_sce))`) as described in {ref}`Additional SingleCellExperiment components for CITE-seq libraries (with ADT tags) <sce_file_contents:Additional SingleCellExperiment components for CITE-seq libraries (with ADT tags)>`.

For `AnnData` objects, these same quality-control statistics are provided in the `obs` slot of the `AnnData` object as described in {ref}`Additional AnnData components for CITE-seq libraries (with ADT tags) <sce_file_contents: Additional AnnData components for CITE-seq libraries (with ADT tags)>`.

We recommend filtering out these low-quality cells before proceeding with downstream analyses.

Here are some additional resources that can be used for working with ADT counts from CITE-seq experiments:
- [Integrating with Protein Abundance, Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.15/OSCA.advanced/integrating-with-protein-abundance.html)
- [Seurat vignette on using with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette.html)


## Special considerations for multiplexed samples

If the dataset that you have downloaded contains samples that were multiplexed (i.e. cells from multiple samples have been combined together into one library), then the steps to take after downloading will be a little different than for libraries that only contain cells corresponding to a single sample.
Here, multiplexed samples refer to samples that have been combined together into one library using cell hashing and then sequenced together.
This means that a single library contains cells or nuclei that correspond to multiple samples.
Each sample has been tagged with a hashtag oligo (HTO) prior to mixing, and that HTO can be used to identify which cells belong to which sample within a multiplexed library.
The libraries available for download on the portal have not been separated by sample (i.e. demultiplexed), and therefore contain data from multiple samples.

Note that multiplexed sample libraries are only available as `SingleCellExperiment` objects, and are not currently available as `AnnData` objects.
If you prefer to work with `AnnData` objects, we recommend using the [`zellkonverter` package](https://bioconductor.org/packages/3.17/bioc/html/zellkonverter.html) to convert the `SingleCellExperiment` object to a HDF5 file containing an `AnnData` object.


Libraries containing multiplexed samples can be initially processed using the same workflow described above including removal of [low quality cells](#quality-control), [normalization](#normalization), and [dimensionality reduction](#dimensionality-reduction).
Demultiplexing can then be used to identify the sample that each cell is from.
Demultiplexing has already been performed using both [`Seurat::HTODemux`](https://satijalab.org/seurat/reference/htodemux) and [`DropletUtils::hashedDrops`](https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html).
For samples where corresponding bulk RNA-sequencing data is available, {ref}`genetic demultiplexing <processing_information:Genetic demultiplexing>` was also conducted.
The associated demultiplexing results were summarized and are available in the `colData` of the processed `SingleCellExperiment` object.
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
