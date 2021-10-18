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
[He _et al._ (2021)](https://doi.org/10.1101/2021.06.29.450377)) demonstrated that Alevin-fry can process single-cell and single-nuclei data more quickly and efficiently then other available methods, while also decreasing the false positive rate of gene detection that is commonly seen in methods that utilize transcriptome alignment.
[You _et al._ (2021)](https://doi.org/10.1101/2021.06.17.448895) and [Tian _et al._ (2019)](https://doi.org/10.1038/s41592-019-0425-8) have also noted that results from different pre-processing workflows for single-cell RNA-sequencing analysis tend to result in compatible results downstream.

## How do I use the provided RDS files in R? 

We are providing the gene expression data to you as a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) in an RDS file.

_Note: You will need to install and load the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor to work with the provided files._

To read in the RDS files you can use the `readRDS` command in base R. 

```
library(SingleCellExperiment)
scpca_sample <- readRDS("SCPCL000000_filtered.rds")
```

## What is the difference between participants, samples, and libraries?

The `participant_id` indicates the participant from which a sample or collection of samples was obtained. 
A sample ID, labeled as `scpca_sample_id` and indicated by the prefix `SCPCS`, represents a unique tissue that was collected from a participant. 
For example, one participant may have a sample collected both at initial diagnosis and at relapse.
This would result in two different sample ID's with the same participant ID. 

The library ID, labeled as `scpca_library_id` and indicated by the prefix `SCPCL`, represents a single set of cells from a tissue sample.
For single-cell or single-nuclei experiments, this will be the result of emulsion and droplet generation using the 10X Genomics workflow, potentially including both RNA-seq, CITE-seq and cell hashing sequencing libraries. 
For a bulk RNA-seq experiment, this will result in a single sequencing library. 

In most cases, each sample will only have one corresponding single-cell library, and may also have an associated bulk RNA-seq library.
However, in some cases multiple libraries underwent separate droplet generation and sequencing from the same sample, resulting in more than one library ID being associated with the same sample ID. 

## What genes are included in the reference transcriptome? 

The {ref}`reference transcriptome index <processing_information:reference transcriptome index>` that was used for alignment was constructed by extracting both spliced cDNA and intronic regions from the primary genome assembly GRCh38, Ensembl database version 104 ([see the code used to generate the reference transcriptome](https://github.com/AlexsLemonade/scpca-nf/blob/main/bin/make_splici_fasta.R)).
The resulting reference transcriptome index contains 60,319 genes.
In addition to protein-coding genes, this list of genes includes pseudogenes and non-coding RNA.
The gene expression data files available for download report all possible genes present in the reference transcriptome, even if not detected in a given library. 

## Where can I see the code for generating QC reports? 

A QC report for every processed library is included with all downloads, generated from the unfiltered and {ref}`filtered <processing_information:filtering cells>` {ref}`gene expression files <gene_expression_file_contents:gene expression file contents>`.
You can find the [function for generating a QC report](https://github.com/AlexsLemonade/scpcaTools/blob/main/R/generate_qc_report.R) and the [QC report template documents](https://github.com/AlexsLemonade/scpcaTools/tree/main/inst/rmd) in the package we developed for working with processed ScPCA data, [`scpcaTools`](https://github.com/AlexsLemonade/scpcaTools). 

## What if I want to use Seurat instead of Bioconductor? 

The files available for download contain [`SingleCellExperiment` objects](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html). 
If desired, these can be converted into Seurat objects. 

You will need to [install and load the `Seurat` package](https://satijalab.org/seurat/articles/install.html) to work with Seurat objects.

For libraries that only contain RNA-sequencing data (i.e. do not have a CITE-seq library found in the `altExp` of the `SingleCellExperiment` object), you can use the following commands:

```
library(Seurat)

# grab counts matrix from the SingleCellExperiment
counts <- as.matrix(counts(sce))

# create seurat object 
seurat_object <- CreateSeuratObject(counts = counts, assay = "RNA",
                                    project = "project_name",
                                    # add colData from SingleCellExperiment
                                    meta.data = as.data.frame(colData(sce)))

# add rowData from SingleCellExperiment to Seurat 
seurat_object[["RNA"]]@meta.features <- as.data.frame(rowData(sce))

# add metadata from SingleCellExperiment to Seurat
seurat_object@misc <- metadata(sce)
```

In order to convert `SingleCellExperiment` objects from libraries that have counts from both RNA-sequencing and CITE-seq, you can use the following additional commands.
This adds a second assay containing the CITE-seq counts and associated feature data to an existing Seurat object:

```
# create assay object in Seurat from CITE-seq counts found in altExp(SingleCellExperiment)
cite_counts <- counts(altExp(sce))
cite_assay <- CreateAssayObject(counts = cite_counts)

# add rowData from altExp to assay 
cite_assay@meta.features = as.data.frame(rowData(altExp(sce)))

# add altExp from SingleCellExperiment as second assay to Seurat
seurat_object[["CITE"]] <- cite_assay
```

## What if I want to use Python instead of R? 

Gene expression data files are only available as RDS files, and in order to view the contents must be opened using R. 
If you prefer to work in python, count data can first be converted into the 10X format using [`DropletUtils::write10xCounts()`](https://rdrr.io/bioc/DropletUtils/man/write10xCounts.html). 
Note that you will need to install the [`DropletUtils` package](https://www.bioconductor.org/packages/devel/bioc/html/DropletUtils.html) to use this function.

This will result in three files being produced:
- the counts matrix in sparse matrix format - `matrix.mtx.gz`
- the row names, or gene names, saved as a TSV - `features.tsv.gz`
- the column names, or cell barcodes, saved as a TSV - `barcodes.tsv.gz`

```
library(SingleCellExperiment)

# read in the RDS file to be converted
sce <- readRDS("SCPCL000000_filtered.rds")

# extract the counts matrix to be saved 
rna_counts <- counts(sce)

# write counts to 10X format
DropletUtils::write10xCounts("SCPCL000000-rna", rna_counts, 
                             barcodes = colnames(rna_counts),
                             gene.id = rownames(rna_counts))

```

If a library has associated CITE-seq that exists, you will have to save that separately.

```
# extract the CITE counts matrix to be saved 
cite_counts <- counts(altExp(sce))

# write counts to 10X format
DropletUtils::write10xCounts("SCPCL000000-cite", cite_counts, 
                             barcodes = colnames(cite_counts),
                             gene.id = rownames(cite_counts))

```

These files can then be directly read into python using the [`scanpy` package](https://scanpy.readthedocs.io/en/stable/), creating an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html).
Note that you will need to [install the `scanpy` package](https://scanpy.readthedocs.io/en/stable/installation.html).

```
import scanpy as sc

#read in 10X formatted files
rna_file_directory = "SCPCL000000-rna"
anndata_object = sc.read_10x_mtx(rna_file_directory)

cite_file_directory = "SCPCL000000-cite"
cite_anndata = sc.read_10x_mtx(cite_file_directory)

# append CITE-seq anndata to RNA-seq anndata
anndata_object['CITE'] = cite_anndata.to_df()
```

It should be noted that in this conversion the `colData`, `rowData` and metadata that is found in the original `SingleCellExperiment` objects will not be obtained. 
If you would like to maintain this data, we provide a [reference from the authors of scanpy](https://theislab.github.io/scanpy-in-R/#converting-from-r-to-python) discussing one method of saving pieces of the `SingleCellExperiment` object separately and adding them into an `AnnData` object.
