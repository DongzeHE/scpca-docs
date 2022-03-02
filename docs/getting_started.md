# Getting started with a ScPCA dataset 

This section provides information on next steps to take after downloading samples or an entire project from the ScPCA portal. 

## Importing ScPCA data into R 

Quantified single-cell or single-nuclei gene expression data is provided as an RDS file as described in the {ref}`single cell gene expression file contents section<sce_file_contents:single-cell gene expression file contents>`. 
There are two RDS files that will be available for each library, a `filtered.rds` and an `unfiltered.rds` file. 
The `filtered.rds` files will contain the gene expression data for only droplets that are considered to contain cells, removing any potential empty droplets. 
The `unfiltered.rds` file will contain the gene expression data for all droplets, regardless of the presence of a cell or not. 
See the {ref}`section on filtering cells<processing_information:filtering cells>` for more information on how we remove potential empty droplets. 
In most scenarios, we recommend starting with the `filtered.rds` file. 

The first step in analyzing the provided gene expression data will be to import the data into R. 
As a reminder, each RDS file contains a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html). 
Refer to {ref}`single-cell gene expression file contents<sce_file_contents:single-cell gene expression file contents>` for a more detailed description on the contents of the included `SingleCellExperiment` object. 

In order to work with `SingleCellExperiment` objects in R, we need to ensure that we have the [`SingleCellExperiment` package]() installed and loaded. 

```r
# if SingleCellExperiment is not installed, install the package
# otherwise the installation step can be skipped
BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)
# read in the RDS file, including the complete path to the file's location
sce <- readRDS("SCPCL000000_filtered.rds")
```

## Filtering my data 

## Downstream analysis

## What if I want to use Seurat? 

