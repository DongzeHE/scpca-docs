# Getting started with a ScPCA dataset 

This section provides information on next steps to take after downloading a dataset from the ScPCA portal. 

## Importing ScPCA data into R 

Quantified single-cell or single-nuclei gene expression data is provided as an RDS file as described in the {ref}`single cell gene expression file contents section<sce_file_contents:single-cell gene expression file contents>`. 
There are two RDS files that will be available for each library, a `filtered.rds` and an `unfiltered.rds` file. 
The `filtered.rds` files will contain the gene expression data for only droplets that are considered to contain cells, removing any potential empty droplets. 
The `unfiltered.rds` file will contain the gene expression data for all droplets, regardless of the presence of a cell or not. 
See the {ref}`section on filtering cells<processing_information:filtering cells>` for more information on how we remove potential empty droplets. 
In most scenarios, we recommend starting with the `filtered.rds` file. 

The first step in analyzing the provided gene expression data will be to import the data into R. 
As a reminder, each RDS file contains a [`SingleCellExperiment` object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html). 
Refer to {ref}`single-cell gene expression file contents<sce_file_contents:single-cell gene expression file contents>` for a more detailed description on the contents of the included `SingleCellExperiment` object. 

More resources for learning about `SingleCellExperiment` objects: 

- [Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) on working with `SingleCellExperiment` objects
- [Orchestrating Single Cell Analysis chapter on the `SingleCellExperiment` class](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html)

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

## Filtering my data 

After we have imported the RDS file into R and accessed the `SingleCellExperiment` object, we can begin working with the data. 
Before we perform any downstream steps, it is recommended to remove any low quality cells from our dataset.
This would include cells that may be dying or dead, showing a higher percentage of reads coming from mitochondrial cells, or cells that may not have been sequenced as thoroughly as others in the same library, showing a lower number of total reads and unique genes identified. 

### Filtering with miQC 

For all filtered `SingleCellExperiment` objects present in the `filtered.rds` files, we have used [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html), a data driven approach to predict low-quality cells, to identify cells that should be removed prior to downstream analyses. 
`miQC` jointly models the proportion of mitochondrial reads and the number of unique genes detected in each cell to calculate the probability of a cell being compromised (i.e. dead or damaged). 
High-quality cells are those with a low probability of being being compromised (< 0.75) or sufficiently low mitochondrial content. 
We have already incorporated the `miQC` results into the {ref}`colData() of the SingleCellExperiment object<sce_file_contents:cell metrics>`. 
All cells that are identified as low-quality cells will have a `FALSE` in the `miQC_pass` column of the `colData()`. 

The following commands can be used to remove those cells prior to downstream analysis: 

```r
# create a vector of barcode names for cells that are not compromised 
# the `$` notation denotes access to the colData slot of the `SingleCellExperiment` object 
cells_to_keep <- colnames(sce$miQC_pass == TRUE)

# filter the `SingleCellExperiment`
filtered_sce <- sce[, cells_to_keep]
```
For more information on using `miQC` for filtering cells, see the [`miQC` vignette.](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html)

### Additional Resources 



## Downstream analysis

## What if I want to use Seurat? 

