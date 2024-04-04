# Frequently Asked Questions

## Why did we use Alevin-fry for processing?

We aimed to process all of the data in the portal such that it is comparable to widely used pipelines, namely Cell Ranger from 10x Genomics.
In our own benchmarking, we found that [Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) produces very similar results to [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count), while allowing faster, more memory efficient processing of single-cell and single-nuclei RNA-sequencing data.
In the configuration that we are using ("selective alignment" mapping to a human transcriptome that includes introns), Alevin-fry uses approximately 12-16 GB of memory per sample and completes mapping and quantification in less than an hour.
By contrast, Cell Ranger uses up to 25-30 GB of memory per sample and takes anywhere from 2-8 hours to align and quantify one sample.
Quantification of samples processed with both Alevin-fry and Cell Ranger resulted in similar distributions of mapped UMI count per cell and genes detected per cell for both tools.

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/total_umi_per_cell.png?raw=true)

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/total_genes_per_cell.png?raw=true)

We also compared the mean gene expression reported for each gene by both methods and observed a high correlation with a Pearson R correlation coefficient of 0.98.

![](https://github.com/AlexsLemonade/alsf-scpca/blob/c0c2442d7242f6e06a5ac6d1e45bd1951780da14/analysis/docs-figures/plots/gene_exp_correlation.png?raw=true)

Recent reports from others support our findings.
[He _et al._ (2021)](https://doi.org/10.1101/2021.06.29.450377) demonstrated that `alevin-fry` can process single-cell and single-nuclei data more quickly and efficiently then other available methods, while also decreasing the false positive rate of gene detection that is commonly seen in methods that utilize transcriptome alignment.
[You _et al._ (2021)](https://doi.org/10.1101/2021.06.17.448895) and [Tian _et al._ (2019)](https://doi.org/10.1038/s41592-019-0425-8) have also noted that results from different pre-processing workflows for single-cell RNA-sequencing analysis tend to result in compatible results downstream.

## How do I use the provided RDS files in R?

If you would like to work with the gene expression data in R, you will need to choose the option for downloading the data as a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
This download includes RDS files that can be directly read into R.

_Note: You will need to install and load the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor to work with the provided files._

To read in the RDS files you can use the `readRDS` command in base R.

```r
library(SingleCellExperiment)
scpca_sample <- readRDS("SCPCL000000_processed.rds")
```

A full description of the contents of the `SingleCellExperiment` object can be found in the section on {ref}`Components of a SingleCellExperiment object <sce_file_contents:Components of a singlecellexperiment object>`.
For more information on working with the RDS files, see {ref}`Getting started with an ScPCA dataset <getting_started:Getting started with an scpca dataset>`.

## How do I use the provided HDF5 files in Python?

If you would like to work with the gene expression data in Python, you will need to choose the option for downloading the data as an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html).
This download includes HDF5 files that can be directly read into Python.

_Note: You will need to install the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html) to work with the provided files._

To read in the HDF5 files you can use the `readh5ad` function from the `AnnData` package.

```python
import anndata
scpca_sample = anndata.readh5ad(file = "SCPCL000000_processed_rna.hdf5")
```

A full description of the contents of the `AnnData` object can be found in the section on {ref}`Components of an AnnData object <sce_file_contents:Components of an anndata object>`.
For more information on working with the HDF5 files, see {ref}`Getting started with an ScPCA dataset <getting_started:Getting started with an scpca dataset>`.

## What is the difference between samples and libraries?

A sample ID, labeled as `scpca_sample_id` and indicated by the prefix `SCPCS`, represents a unique tissue that was collected from a participant.

The library ID, labeled as `scpca_library_id` and indicated by the prefix `SCPCL`, represents a single set of cells from a tissue sample, or a particular combination of samples in the case of multiplexed libraries.
For single-cell or single-nuclei experiments, this will be the result of emulsion and droplet generation using the 10x Genomics workflow, potentially including both RNA-seq, ADT (i.e., from CITE-seq), and cell hashing sequencing libraries.
Multiplexed libraries will have more than one sample ID corresponding to each library ID.

In most cases, each sample will only have one corresponding single-cell or single-nuclei library, and may also have an associated bulk RNA-seq library.
However, in some cases multiple libraries were created by separate droplet generation and sequencing from the same sample, resulting in more than one single-cell or single-nuclei library ID being associated with the same sample ID.

## What is a participant ID?

The `participant_id` is a unique ID provided by the submitter to indicate the participant from which a collection of samples was obtained.
For example, one participant may have a sample collected both at initial diagnosis and at relapse.
This would result in two different sample ID's, but the same participant ID.
However, for most participants, only a single sample was collected and submitted for sequencing.

## What is a multiplexed sample?

Multiplexed samples refer to samples that have been combined together into a single library using cell hashing ([Stoeckius _et al._ 2018](https://doi.org/10.1186/s13059-018-1603-1)) or a related technology and then sequenced together.
This means that a single library contains cells or nuclei that correspond to multiple samples.
Each sample has been tagged with a hashtag oligo (HTO) prior to mixing, and that HTO can be used to identify which cells or nuclei belong to which sample within a multiplexed library.
The libraries available for download on the portal have not been separated by sample (i.e. demultiplexed), and therefore contain data from multiple samples.
For more information on working with multiplexed samples, see the {ref}`special considerations for multiplexed samples section in getting started with an ScPCA dataset <getting_started:Special considerations for multiplexed samples>`.

## Why are demultiplexed samples not available?

Downloading a multiplexed sample on the portal will result in obtaining the gene expression files corresponding to the library containing the chosen multiplexed sample and any other samples that were multiplexed with that chosen sample.
This means that users will receive the gene expression data for all samples that were combined into a given library and will have to separate any cells corresponding to the sample of interest before proceeding with downstream analysis.

We have applied multiple {ref}`demultiplexing methods <processing_information:hto demultiplexing>` to multiplexed libraries and noticed that these demultiplexing methods can vary both in calls and confidence levels assigned.
[Here we have performed some exploratory analysis comparing demultiplexing methods within a single multiplexed library](https://htmlpreview.github.io/?https://github.com/AlexsLemonade/alsf-scpca/blob/main/analysis/quantifier-comparisons/15-demux-comparisons.nb.html).
Because of the inconsistency across demultiplexing methods used, the choice of demultiplexing method to use is at the discretion of the user.
Rather than separating out each sample, the sample calls and any associated statistics regarding sample calls for multiple demultiplexing methods can be found in the `_filtered.rds` file for each multiplexed library.
See the {ref}`demultiplexing results section <sce_file_contents:demultiplexing results>` for instructions on how to access the demultiplexing results in the `SingleCellExperiment` objects for multiplexed libraries.
We also include the Hash Tag Oligo counts matrix to allow demultiplexing using other available methods.

## What are estimated demux cell counts?

Estimated demux cell counts are provided for multiplexed libraries and refer to the estimated cell counts for each sample that is present in the library.
In order to provide an estimate of the number of cells or nuclei that are present in a given sample before download, we use the estimated number of cells per sample identified by one of the tested demultiplexing methods.
However, these estimated demux cell counts should only be considered a guide; we encourage users to investigate the data on their own and make their own decisions on the best demultiplexing method to use for their research purposes.

Note that not all cells in a library are included in the estimated demux cell count, as some cells may not have been assigned to a sample.
Estimated demux cell counts are only reported for multiplexed samples and are not reported for single-cell or single-nuclei samples that are not multiplexed.
For more about demultiplexing, see the section on {ref}`processing multiplexed libraries <processing_information:multiplexed libraries>`.

## What genes are included in the reference transcriptome?

The {ref}`reference transcriptome index <processing_information:reference transcriptome index>` that was used for alignment was constructed by extracting both spliced cDNA and intronic regions from the primary genome assembly GRCh38, Ensembl database version 104 ([see the code used to generate the reference transcriptome](https://github.com/AlexsLemonade/scpca-nf/blob/main/bin/make_reference_fasta.R)).
The resulting reference transcriptome index contains 60,319 genes.
In addition to protein-coding genes, this list of genes includes pseudogenes and non-coding RNA.
The gene expression data files available for download report all possible genes present in the reference transcriptome, even if not detected in a given library.

## Where can I see the code for generating QC reports?

A QC report for every processed library is included with all downloads, generated from the unfiltered and {ref}`filtered <processing_information:filtering cells>` {ref}`Single-cell gene expression files <sce_file_contents:Single-cell gene expression file contents>`.
You can find the [function for generating a QC report](https://github.com/AlexsLemonade/scpcaTools/blob/main/R/generate_qc_report.R) and the [QC report template documents](https://github.com/AlexsLemonade/scpcaTools/tree/main/inst/rmd) in the package we developed for working with processed ScPCA data, [`scpcaTools`](https://github.com/AlexsLemonade/scpcaTools).

## What if I want to use Seurat instead of Bioconductor?

The RDS files available for download contain [`SingleCellExperiment` objects](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
If desired, these can be converted into Seurat objects.

You will need to [install and load the `Seurat` package](https://satijalab.org/seurat/articles/install.html) to work with Seurat objects.

For libraries that only contain RNA-seq data (i.e., do not have an ADT library found in the `altExp` of the `SingleCellExperiment` object), you can use the following commands:

```r
library(Seurat)
library(SingleCellExperiment)

# read in RDS file
sce <- readRDS("SCPCL000000_filtered.rds")

# create seurat object from the SCE counts matrix
seurat_object <- CreateSeuratObject(counts = counts(sce),
                                    assay = "RNA",
                                    project = "SCPCL000000")
```
The above code will only maintain information found in the original counts matrix from the `SingleCellExperiment`.
Optionally, if you would like to keep the included cell and gene associated metadata during conversion to the Seurat object you can perform the below additional steps:

```r
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(colData(sce))
row_metadata <- as.data.frame(rowData(sce))

# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_object@meta.data <- cell_metadata

# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_object[["RNA"]]@meta.features <- row_metadata

# add metadata from SingleCellExperiment to Seurat
seurat_object@misc <- metadata(sce)
```

For `SingleCellExperiment` objects from libraries with both RNA-seq and ADT data, you can use the following additional commands to add a second assay containing the ADT counts and associated feature data:

```r
# create assay object in Seurat from ADT counts found in altExp(SingleCellExperiment)
adt_assay <- CreateAssayObject(counts = counts(altExp(sce)))

# optional: add row metadata (rowData) from altExp to assay
adt_row_metadata <- as.data.frame(rowData(altExp(sce)))
adt_assay@meta.features <- adt_row_metadata

# add altExp from SingleCellExperiment as second assay to Seurat
seurat_object[["ADT"]] <- adt_assay
```

## When should I download a project as a merged object?

When you download all data for a ScPCA project, you will be presented with two options.
You can either download the project such that the data for each sample is stored in separate files, or you can download a single file that contains a merged object with data from all samples in the project.
This merged object contains combined data from all samples (and therefore all libraries), including expression count matrices and associated metadata.
The samples have simply been merged into a single file - _they have not been integrated/batch-corrected_.

You may prefer to download this merged object instead of individual sample files to facilitate downstream analyses that consider multiple samples at once, such as differential expression analysis, integrating multiple samples, or jointly clustering multiple samples.

Please refer to {ref}`the getting started with a merged object section<getting_started:Working with a Merged ScPCA object>` for more details on working with merged objects.


## Which projects can I download as merged objects?

Most projects in the ScPCA Portal are available for download as a merged object.
There are two types of projects for which merged objects are not available:

- Projects comprised of spatial transcriptomics
    - As described in {ref}`the spatial transcriptomics processing section<processing_information:spatial transcriptomics>`, no post-processing is performed on these libraries after running Space Ranger.
    Therefore, merging samples into a single object is beyond the scope of the ScPCA pipeline.

- Projects containing multiplexed libraries
    - Although the ScPCA pipeline {ref}`reports demultiplexing results<processing_information:HTO demultiplexing>`, it does not actually perform demultiplexing.
    As there is no guarantee that a unique HTO was used for each sample in a given project, it would not necessarily be possible to determine which HTO corresponds to which sample in a merged object.

- Projects containing more than 50 samples
    - The more samples that are included in a merged object, the larger the object, and the more difficult it will be to work with that object in R or Python.
    Because of this, we do not provide merged objects for projects with more than 50 samples as the size of the merged object is too large.

## Why doesn't my existing code work on a new download from the Portal?

Although we try to maintain backward compatibility, new features added to the ScPCA Portal may result in downloads that are no longer compatible with code written with older downloads from the ScPCA Portal in mind.
Please see our {ref}`CHANGELOG <CHANGELOG:CHANGELOG>` for a summary of changes that impact downloads from the Portal.
