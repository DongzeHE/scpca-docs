# Frequently Asked Questions 

#### Why did we use Alevin-fry for processing? 

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

#### How do I use the provided RDS files in R? 

We are providing the gene expression data to you as a [`SingleCellExperiment`](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) in an RDS file.

_Note: You will need to install and load the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.13/bioc/html/SingleCellExperiment.html) from Bioconductor to work with the provided files._

To read in the RDS files you can use the following options. 
In base R, you can use the `readRDS` command. 
```
scpca_sample <- readRDS("SCPCL000001_filtered.rds")
```

Or in the `readr` package, you can use the `read_rds()` function. 
```
scpca_sample <- readr::read_rds("SCPCL000001_filtered.rds")
``` 

The `SingleCellExperiment` holds the gene expression data as a matrix in the assays slot and any associated metadata. 

To access the counts data directly you can do: 
```
counts(scpca_sample)
```

We have included metrics for each cell (associated with every column) and each gene (aassociated with each row). 
To access the this data you can do: 
```
# get data associated with each cell barcode
colData(scpca_sample)

# get data associated with each row barcode
rowData(scpca_sample)
```

To access any additional metadata, use: 
```
metadata(scpca_sample)
```