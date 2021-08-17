# Processing Information

For every cDNA library that was submitted to the ScPCA portal, FASTQ files were used as input to our processing pipeline to generate both an unfiltered and filtered gene by cell counts matrix. 
All single-cell, single-nuclei, CITE-seq, and spatial transcriptomics libraries were processed using [Alevin-fry](https://alevin-fry.readthedocs.io/en/latest/)[He _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1.full.pdf).
Using a subset of single-cell and single-nuclei samples, we compared use of Alevin-fry to 10X Genomics' [Cellranger count.](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)
Although Cellranger is more commonly used, we found that Alevin-fry was more computationally efficient, using only 12-16 GB of memory per sample and completing in less than 1 hour. 
See here for a detailed comparison comparison of cell metrics across processing tools revealing minimal differences between Alevin-fry and Cellranger, supporting the use of Alevin-fry in building the ScPCA portal. (include link to a notebook with key benchmarking metrics)

Bulk RNA-sequencing data was aligned and quantified using [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html).
Gene by cell counts matrices were then transformed into [`SingleCellExperiment` objects](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html). 

## Single Cell and Single Nuclei Samples

### Alignment and Quantification using Alevin-fry

FASTQ files from single-cell and single-nuclei RNA-sequencing libraries were first aligned to an index containing transcripts from both spliced and unspliced cDNA, or the `splici` index, using the selective alignment strategy invoked by `salmon alevin`. 
We then obtained the list of all possible cells by using `alevin-fry generate-permit-list --unfiltered-pl`. 
The `--unfiltered-pl` option was used alongside the list of barcodes obtained from 10X Genomics based on the 10X version utilized in the original single-cell or single-nuclei RNA-sequencing library generation. 
Use of this option with alevin-fry returns all possible cells within 1 edit distance from the provided 10X whitelist in the counts matrix. 
Lastly, feature quantification is performed and UMI's are de-duplicated using Alevin-fry's `cr-like-em` resolution strategy. 
The `cr-like-em` resolution strategy assigns all UMI's that align to a single gene to that gene. 
Any gene that is assigned to more than one gene, the UMI is assigned to the gene with the highest count. 
If a UMI is assigned to more than one gene with equal frequency, then the counts for each gene are determined via an expectation maximization algorithm. 
The unfiltered gene by cell counts matrix contains all possible cells found in the 10X whitelist. 

#### The Splici Index

The `splici index` was chosen to be used for single-cell RNA-sequencing libraries because of recent publications highlighting that alignment to the transcriptome alone can lead to spuriously detected genes. 
The rationale for this is that some small fraction of reads from single cell experiments are derived from intronic or intergenic sequnces and not from spliced cDNA [Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full#sec-5). 
To minimize the increase in false positive rates, an [index with both transcripts from spliced cDNA and intronic regions can be used](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) [He _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1.full.pdf). 
In our hands, we have found that use of the `splici` index led to a more comparable distribution of unique genes found per cell to Cellranger than use of an index obtained from spliced cDNA transcripts. 

For single-nuclei RNA-sequencing experiments, we expect reads to be from both nascent mRNA and mature mRNA.
Therefore, we use the `splici` index to allow for quantification of both mature, spliced mRNA and nascent, unspliced mRNA that has been captured. 

##### Generating a Splici Index

In order to generate the `splici` index...

#### Selective Alignment

Alevin-fry utilizes the same selective alignment strategy invoked by `salmon` to align reads. 
Briefly, selective alignment uses a mapping score validated approach to identify maximal exact matches between reads and the provided index. 
Alevin-fry provides the option to use either selective alignment or pseudoalignment, a light weight alignment strategy that does not use mapping scores to validate alignments. 
We observed that use of selective alignment with the `splici` index mitigated any spurious increases in total RNA or total genes per cell found when using alevin-fry with pseudoalignment in comparison to Cellranger.
For all single-cell and single-nuclei samples, selective alignment to the `splici` index was used. 

A more detailed description of the alignment strategy invoked by `salmon` in conjuction with alevin-fry can be found [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8). 

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