# Processing information

## Single-cell and single-nuclei RNA-seq

### Alignment and quantification using alevin-fry

We used [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/) to generate gene by cell counts matrices for all single-cell and single-nuclei samples.
In brief, we utilized [selective alignment](#selective-alignment) to the [`splici` index](#reference-transcriptome-index) for all single-cell and single-nuclei samples. 

#### Reference transcriptome index

For all samples, we aligned FASTQ files to a reference transcriptome index referred to as the `splici` index.
The [`splici` index](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) is built using transcripts from both spliced cDNA and intronic regions.
Inclusion of intronic regions in the index used for alignment allowed us to capture both reads from mature, spliced cDNA and nascent, unspliced cDNA. 
Alignment of RNA-sequencing data to an index containing intronic regions has been shown to reduce spuriously detected genes ([He _et al._ 2021](https://doi.org/10.1101/2021.06.29.450377), [Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full#sec-5))
In our hands, we have found that use of the `splici` index led to a more comparable distribution of unique genes found per cell to Cell Ranger than use of an index obtained from spliced cDNA transcripts only. 

#### Selective alignment

We aligned reads to the transcriptome index using `salmon` with the default "selective alignment" strategy. 
Briefly, selective alignment uses a mapping score validated approach to identify maximal exact matches between reads and the provided index. 
For all samples, we used selective alignment to the `splici` index. 

A more detailed description of the alignment strategy invoked by `salmon` in conjunction with `alevin-fry` can be found in [Srivastava _et al._ 2020](https://doi.org/10.1186/s13059-020-02151-8).

#### Alevin-fry parameters 

After we aligned FASTQ files using selective alignment to the `splici` index, we continued with the `alevin-fry` pipeline using the following parameters: 

1. During the [`generate-permit-list` step of `alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/generate_permit_list.html), we used the `--unfiltered-pl` option, which returns any cell with at least 1 read found in a reference barcode list. 
For our reference barcode list, we used a list of all possible cell barcodes from 10X Genomics.

2. We chose to use the `cr-like-em` resolution strategy for [feature quantification and UMI de-duplication](https://alevin-fry.readthedocs.io/en/latest/quant.html). 
Similar to the way Cell Ranger performs feature quantification, the `cr-like-em` resolution strategy assigns all UMIs that align to a single gene to that gene. 
In contrast to Cell Ranger, `cr-like-em` keeps multi-mapped reads and invokes an extra step to assign these multi-mapped reads to a UMI.

### Post alevin-fry processing

#### Combining counts from spliced cDNA and intronic regions

For single-cell samples, we only included reads aligning to spliced cDNA transcripts in the counts matrix. 
For single-nuclei samples, all counts for spliced cDNA and intronic regions were summed for each gene to return the total counts summarized by gene in the counts matrix. 

#### Filtering cells

In addition to an unfiltered counts matrix, we provide a matrix filtered to only droplets that are likely to include cells, removing any cell barcodes that are likely to correspond to empty droplets.
To do this we used [`DropletUtils::emptyDrops()`](https://www.bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets), a function that estimates the profile of cells containing ambient RNA and tests the likelihood of all other droplets as differing from the ambient profile [Lun _et al._ 2019](https://doi.org/10.1186/s13059-019-1662-y). 
We defined the ambient profile to include all droplets with less than 200 UMI per cell by using the option `lower=200`.
We consider droplets with an FDR less than or equal to 0.01 to be cell-containing droplets. 
Only cells that pass this FDR threshold are included in the filtered counts matrix.

### CITE-seq quantification

#### Alignment and quantification using alevin-fry

#### Combining CITE counts with RNA counts