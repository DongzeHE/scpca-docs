# Getting started with a ScPCA dataset 

This section provides information on next steps to take after downloading a dataset from the ScPCA portal. 

## Importing ScPCA data into R 

Quantified single-cell or single-nuclei gene expression data is provided as an RDS file as described in the {ref}`single cell gene expression file contents section<sce_file_contents:single-cell gene expression file contents>`. 
There are two RDS files that will be available for each library, a `filtered.rds` and an `unfiltered.rds` file. 
The `filtered.rds` files will contain the gene expression data for only droplets that are likely to contain cells, removing any potential empty droplets. 
The `unfiltered.rds` file will contain the gene expression data for all droplets, regardless of the presence of a cell or not. 
See the {ref}`section on filtering cells<processing_information:filtering cells>` for more information on how we remove potential empty droplets. 
In most scenarios, we recommend starting with the `filtered.rds` file. 

The first step in analyzing the provided gene expression data will be to import the data into R. 
As a reminder, each RDS file contains a [`SingleCellExperiment` object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html). 
Refer to {ref}`single-cell gene expression file contents<sce_file_contents:single-cell gene expression file contents>` for a more detailed description on the contents of the included `SingleCellExperiment` object. 

In order to work with `SingleCellExperiment` objects in R, we need to ensure that we have the [`SingleCellExperiment` package]() installed and loaded. 

The below commands can be used to import the RDS file into R and save the `SingleCellExperiment` object. 

```r
# if SingleCellExperiment is not installed, install the package
# otherwise the installation step can be skipped
BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)
# read in the RDS file, including the complete path to the file's location
sce <- readRDS("SCPCL000000_filtered.rds")
```

More resources for learning about `SingleCellExperiment` objects: 

- [`SingleCellExperiment` objects Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) on working with `SingleCellExperiment` objects
- [Orchestrating Single Cell Analysis chapter on the `SingleCellExperiment` class](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html)

## Quality control  

After we have imported the RDS file into R and accessed the `SingleCellExperiment` object, we can begin working with the data. 
Before we perform any downstream steps, it is recommended to remove any low quality cells from our dataset.
This would include cells that may be dying or damaged, showing a higher percentage of reads coming from mitochondrial cells. 
Low quality cells may also be present due to technical artifiacts such as inefficient reverse transcription or PCR amplification resulting in a lower number of total reads and unique genes identified. 

For the filtered `SingleCellExperiment` object present in all `filtered.rds` files, we have used [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html), a data driven approach to predict low-quality cells, to identify cells that should be removed prior to downstream analyses. 
`miQC` jointly models the proportion of mitochondrial reads and the number of unique genes detected in each cell to calculate the probability of a cell being compromised (i.e. dead or damaged). 
High-quality cells are those with a low probability of being being compromised (< 0.75) or sufficiently low mitochondrial content. 
We have already incorporated the `miQC` results into the {ref}`colData() of the SingleCellExperiment object<sce_file_contents:cell metrics>`. 
All cells that are identified as low-quality cells will have `FALSE` in the `miQC_pass` column of the `colData()`. 

The following commands can be used to remove those cells prior to downstream analysis: 

```r
# create a vector of barcode names for cells that are not compromised 
# the `$` notation denotes access to the colData slot of the `SingleCellExperiment` object 
cells_to_keep <- colnames(sce$miQC_pass == TRUE)

# filter the `SingleCellExperiment`
filtered_sce <- sce[, cells_to_keep]
```
For more information on using `miQC` for filtering cells, see the following resources: 
- [`miQC` vignette](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html)
- [Hippen _et al._ 2021](https://doi.org/10.1371/journal.pcbi.1009290)

Alternatively, you can directly filter cells based on having a minimum number of unique genes, total reads, and a maximum fraction of mitochondrial reads. 
Using the function, [`scuttle::addPerCellQCMetrics()`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html), we can calculate the total RNA counts, unique genes, and mitochondrial fraction for each cell and store the results in the `colData()`. 
These metrics can be used to directly filter the `SingleCellExperiment` object based on informed minimum thresholds. 
The `filtered.rds` files contain `SingleCellExperiment` objects with these metrics added to the `colData()` and you can read more about all metrics that are included in the {ref}`cell metrics section of the SingleCellExperiment file contents<sce_file_contents:cell metrics>`.

If you are planning to filter low quality cells using these hard cutoffs, we encourage you to read more about the various metrics and plot the distribution of each metric before deciding on which cells to exclude. 
The chapter on [Quality Control in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/quality-control.html#quality-control) provides a nice guide to checking diagnostic plots and then choosing cutoffs. 

## Normalization 

Single-cell RNA sequencing inherently comes with technical confounders that, if not accounted for, can impact the interpretation of downstream analysis. 
Technical variables include differences in sequencing depth and cell size across cells in a single-cell library ([Stegle _et al._ 2015](https://doi.org/10.1038/nrg3833).)

If we were to proceed with downstream analysis such as differential gene expression and cell clustering without accounting for such effects, these technical problems could impact the overall conclusions.
We need to normalize our data to account for the technical variation between cells in our data and maximize biological variation, allowing us to make more accurate conclusions from downstream analysis.

We recommend using the `scran` and `scater` packages to estimate size factors based on a pool or cluster of cells ([Lun _et al._ 2016](https://doi.org/10.1186/s13059-016-0947-7)).
Here, before calculating size factors, cells are first pooled into clusters based on similar patterns of gene expression using an approximation of principal component analysis (PCA). 
Once cells are assigned to a cluster, size factors are computed for each pool first and then "deconvolved" across each cell in the pool. 
You can read more about implementing this method from the [Normalization by deconvolution chapter in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/normalization.html#normalization-by-deconvolution). 

Follow these steps to perform normalization using `scran` and `scater`: 

```r
# Cluster similar cells 
qclust <- scran::quickCluster(filtered_sce)

# Compute sum factors for each cell cluster grouping.  
filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)

# Normalize and log transform. 
normalized_sce <- scater::logNormCounts(filtered_sce)
```

Here we provide more resources on understanding Normalization in single-cell RNA-seq analysis: 

- [Orchestrating Single Cell Analysis Chapter on Normalization](http://bioconductor.org/books/3.14/OSCA.basic/normalization.html)
- [Hemberg lab scRNA-seq course section on Normalization methods](https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html#normalization-theory)
- Review on Computational challenges in single-cell, including a summary on Normalization and technical variance in scRNA-seq([Stegle _et al._ 2015](https://doi.org/10.1038/nrg3833)).

## Downstream analysis

### Feature selection

Now that our data is normalized we are ready to perform downstream analysis such as dimensionality reduction and clustering. 
In performing both of these procedures, we want to be sure that we are maximizing biological variance and decreasing technical noise. 
To do this we want to select the most highly variable genes and use those for dimensionality reduction. 

```r
# calculate the gene variance 
gene_variance <- scran::modelGeneVar(normalized_sce)
# select the most variable genes
highvar_genes <- scran::getTopHVGs(gene_variance, n = 2000)
```

Visit the chapter on [Feature Selection in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/feature-selection.html#feature-selection) to read more about modeling gene variance and selecting the highly variable genes. 

### Dimensionality Reduction 

In single-cell RNA seq, every gene is another dimension and visualizing high dimension data is problematic. 
It is assumed that genes that are affected by the same biological process are correlated, therefore we do not need to store each gene as an individual dimension. 
Dimensionality reduction is commonly used to reduce dimensions used for plotting, clustering, and other downstream analysis. 

We recommend first performing [principal component analysis (PCA)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#principal-components-analysis), a technique that identifies new axes that capture the largest amount of variation in the data. 
The PCA results can be calculated and stored in the `SingleCellExperiment` object using the following command: 

```r
# calculate a PCA matrix using the most highly variable genes
normalized_sce <- runPCA(normalized_sce, subset_row = highvar_genes)
```

PCA is a linear dimensionality reduction technique and therefore it cannot efficiently pack differences in d dimensions into the first 2 principal components which are used for visualizations. 
Therefore we need a non linear method, such as [UMAP (Uniform Manifold Approximation and Projection)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#uniform-manifold-approximation-and-projection), for visualization. 
UMAP allows for better separation between clusters of cells, but can be dependent on the choice of parameters, such as the number of neighbors and minimum distance between points. 
It's important to note that while the observed clusters do have some meaning, the distance between clusters and the cluster density usually is not related to the similarity or dissimilarity of the clusters. 
Additionally, if the results are completely dependent on the choice of parameters then you should interpret the results with caution. 

UMAP can also be quite slow for a large dataset, so it's recommended to run UMAP with the top 50 principal components, speeding up the analysis. 

```r
# Run UMAP using already stored PCA results
normalized_sce <- runUMAP(normalized_sce, 
                          dimred = "PCA")
```

See below for more resources on dimensionality reduction: 

- [Dimensionality Reduction in Orchestrating Single Cell Analysis](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#dimensionality-reduction)
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
- [HBC Course on single-cell RNA-seq analysis with Seurat](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html) 
