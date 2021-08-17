# Processing Information

## Single Cell and Single Nuclei Samples

### Alignment and Quantification using Alevin-fry

We used [Alevin-fry](https://alevin-fry.readthedocs.io/en/latest/) to generate gene by cell counts matrices for all single-cell and single-nuclei samples.
In brief, we utilized [selective alignment](####selective-alignment) to the [`splici` index](####the-splici-index) for all single-cell and single-nuclei samples. 

#### The Splici Index

The [`splici` index](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) is built using transcripts from both spliced cDNA and intronic regions.
Inclusion of intronic regions in the index used for alignment, allowed us to capture both reads from mature, spliced cDNA and nascent, unspliced cDNA, commonly observed in RNA-sequencing generated from single-nuclei. 
Although we only anticipate nascent, unspliced cDNA to be present in single-nuclei samples, we chose to align the FASTQ files from both single-cell and single-nuclei samples to the `splici` index.
Recent publications highlight that a small fraction of reads from single cell experiments are derived from intronic or intergenic sequences and not from spliced cDNA [Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full#sec-5).
Alignment to the transcriptome alone, without inclusion of intronic regions, can result in an inaccurate assignment of the intronic sequences to spliced cDNA that share sequence similarity with a spliced transcript. 
This has been observed to lead to an increase in spuriously detected genes [He _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1.full.pdf).
In our hands, we have found that use of the `splici` index led to a more comparable distribution of unique genes found per cell to Cell Ranger than use of an index obtained from spliced cDNA transcripts. 

#### Selective Alignment

Alevin-fry utilizes the selective alignment strategy invoked by `salmon` to align reads. 
Briefly, selective alignment uses a mapping score validated approach to identify maximal exact matches between reads and the provided index. 
For all samples, we used selective alignment to the `splici` index. 

A more detailed description of the alignment strategy invoked by `salmon` in conjuction with alevin-fry can be found [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8). 

#### Alevin-Fry Parameters 

After we aligned FASTQ files using selective alignment to the `splici` index, we continued with the alevin-fry pipeline using the following parameters: 

1. During the `generate-permit-list` step of Alevin-fry, we utilized the `--unfiltered-pl` option which returns any cell with at least 1 read found in the reference barcode list that is used. 
For our reference barcode list, we used a list of all possible cell barcodes from 10X Genomics.
Use of this option allowed us to obtain an unfiltered counts matrix with all possible cell barcodes detected in the experiment. 

2. Alevin-fry provides multiple options for performing feature quantification and UMI de-duplication. 
We chose to use the `cr-like-em` resolution strategy for UMI de-duplication. 
The `cr-like-em` resolution strategy assigns all UMI's that align to a single gene to that gene. 
Any gene that is assigned to more than one gene, the UMI is assigned to the gene with the highest count. 


#### Collapsing Reads from spliced cDNA and intronic regions

In building the `splici` index to be used with both single-cell and single-nuclei samples, the genes were split into unspliced (or intronic) regions of the gene and spliced regions of the gene. 
Gene names were then tagged with a `-U` representing unspliced regions or `-S` for spliced regions in the counts matrix output from alevin-fry. 
If a sample was single-cell RNA-sequencing, any counts aligning to gene names ending in `-U` were discarded before returning the final counts matrix. 
If a sample was single-nuclei RNA-sequencing, all counts were retained and counts for intronic and spliced regions of genes were combined to obtain the total for that gene and return the final counts matrix. 

### Filtering Cells

### Samples with CITE-seq

#### Alignment and Quantification using Alevin-fry

#### Combining CITE counts with RNA counts

### Samples with Spatial Transcriptomics

#### Alignment and Quantification using Alevin-fry

## Bulk RNA Samples

### Alignment and Quantification using Salmon

### Normalization of RNA counts