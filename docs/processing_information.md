# Processing information

## Single-cell and single-nuclei RNA-seq

### Mapping and quantification using alevin-fry

We used [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/) to generate gene by cell counts matrices for all single-cell and single-nuclei samples.
In brief, we utilized [selective alignment](#selective-alignment) to the [`splici` index](#reference-transcriptome-index) for all single-cell and single-nuclei samples. 

#### Reference transcriptome index

For all samples, we aligned FASTQ files to a reference transcriptome index referred to as the `splici` index.
The [`splici` index](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) is built using transcripts from both spliced cDNA and intronic regions.
Inclusion of intronic regions in the index used for alignment allowed us to capture both reads from mature, spliced cDNA and nascent, unspliced cDNA. 
Alignment of RNA-seq data to an index containing intronic regions has been shown to reduce spuriously detected genes ([He _et al._ 2021](https://doi.org/10.1101/2021.06.29.450377), [Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full#sec-5))
In our hands, we have found that use of the `splici` index led to a more comparable distribution of unique genes found per cell to Cell Ranger than use of an index obtained from spliced cDNA transcripts only. 

#### Selective alignment

We mapped reads to the transcriptome index using `salmon` with the default "selective alignment" strategy. 
Briefly, selective alignment uses a mapping score validated approach to identify maximal exact matches between reads and the provided index. 
For all samples, we used selective alignment to the `splici` index. 

A more detailed description of the mapping strategy invoked by `salmon` in conjunction with `alevin-fry` can be found in [Srivastava _et al._ 2020](https://doi.org/10.1186/s13059-020-02151-8).

#### Alevin-fry parameters 

After mapping FASTQ files using selective alignment to the `splici` index, we continued with the `alevin-fry` pipeline using the following parameters: 

1. During the [`generate-permit-list` step of `alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/generate_permit_list.html), we used the `--unfiltered-pl` option, which returns any cell with at least 1 read found in a reference barcode list. 
For our reference barcode list, we used a list of all possible cell barcodes from 10X Genomics.

2. We chose to use the `cr-like-em` resolution strategy for [feature quantification and UMI de-duplication](https://alevin-fry.readthedocs.io/en/latest/quant.html). 
Similar to the way Cell Ranger performs feature quantification, the `cr-like-em` resolution strategy assigns all UMIs that align to a single gene to that gene.
In contrast to Cell Ranger, `cr-like-em` keeps multi-mapped reads and invokes an extra step to assign these multi-mapped reads to a UMI.

3. With initial mapping to the `splici` index, `alevin-fry` quantification resulted in separate counts for spliced and unspliced transcripts, and an ambiguous count for reads compatible with either spliced or unspliced transcripts.

### Post alevin-fry processing

#### Combining counts from spliced cDNA and intronic regions

For single-cell samples, we only included reads aligning to spliced and ambiguous cDNA transcripts in the counts matrix. 
For single-nuclei samples, all counts for spliced cDNA and intronic regions were summed for each gene to return the total counts summarized by gene in the counts matrix. 

After combining read counts, values are rounded to integer values.

#### Filtering cells

In addition to an unfiltered counts matrix, we provide a matrix filtered to only cell barcodes from droplets that are likely to include true cells.
To do this we used [`DropletUtils::emptyDropsCellRanger()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDropsCellRanger.html), a function that estimates the profile of cells containing ambient RNA and tests the likelihood of all other droplets as differing from the ambient profile [Lun _et al._ 2019](https://doi.org/10.1186/s13059-019-1662-y). 
This function more closely mimics the filtering performed in Cell Ranger than its predecessor [`DropletUtils::emptyDrops()`](https://www.bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets). 
We consider droplets with an FDR less than or equal to 0.01 to be cell-containing droplets. 
Only cells that pass this FDR threshold are included in the filtered counts matrix.

For some libraries, `DropletUtils::emptyDropsCellRanger()` may fail due to low numbers of droplets with reads or other violations of its assumptions.
For these libraries, only droplets containing at least 100 UMI are included in the filtered counts matrix.

## CITE-seq quantification

CITE-seq libraries with reads from antibody-derived tags (ADTs) were also quantified using  [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/), rounded to integer values.

Reference indices were constructed from the submitter-provided list of antibody barcode sequences corresponding to each library using the `--features` flag of `salmon index`.
Mapping to these indices followed the same procedures as for RNA-seq data, including mapping with [selective alignment](#selective-alignment) and subsequent [quantification via alevin-fry](#alevin-fry-parameters).

### Combining CITE counts with RNA counts

The unfiltered CITE-seq and RNA-seq count matrices often include somewhat different sets of cell barcodes, due to stochastic variation in library construction and sequencing.
When normalizing these two count matrices to the same set of cells, we chose to prioritize RNA-seq results for broad comparability among libraries with and without CITE-seq data. 
Any cell barcodes that appeared only in CITE-seq data were discarded.
Cell barcodes that were present only in the RNA-seq data (i.e., did _not_ appear in the CITE-seq data) were assigned zero counts for all ADTs. 
When cells were [filtered based on RNA-seq content](#filtering-cells) after quantification, the CITE-seq count matrix was filtered to match.

## Spatial Transcriptomics

**_Data Available Soon_**

### Mapping and quantification using Space Ranger

Processing spatial transcriptomics libraries requires two steps - gene expression quantification and tissue detection. 
In the absence of independent tissue detection methods to use with Alevin-fry, we used [10X Genomics' Space Ranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger) to obtain both gene expression and spatial information. 
[`spaceranger count`](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count) takes FASTQ files and a microscopic slide image as input and performs alignment, quantification, and tissue detection for each spot. 
In contrast to Alevin-fry, which maps reads to a [reference transcriptome index](#reference-transcriptome-index), Space Ranger aligns transcript reads to the reference genome using STAR ([Dobin _et al_., 2012](https://doi.org/10.1093/bioinformatics/bts635)).
See the 10X documentation for more information on how Space Ranger [quantifies gene expression](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/algorithms/overview) and [detects tissue](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/algorithms/imaging).

## Bulk RNA samples 

### Preprocessing with fastp

Prior to quantifying gene expression for bulk RNA-seq samples, FASTQ files are pre-processed using [`fastp`](https://github.com/OpenGene/fastp) to perform adapter trimming, quality filtering, and length filtering. 
For length filtering, trimmed reads shorter than 20 basepairs were removed by using the `--length_required 20` option. 
All other filtering and trimming was performed using the default strategies enabled in `fastp`.  

### Mapping and quantification using salmon

To quantify gene expression for bulk RNA-seq samples, we used [`salmon quant`](https://salmon.readthedocs.io/en/latest/salmon.html).
Here, we performed selective alignment of the trimmed and filtered FASTQ files to a decoy-aware reference transcriptome index ([Srivastava _et al._ 2020](https://doi.org/10.1186/s13059-020-02151-8)). 
The [decoy-aware reference transcriptome](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode), was created from spliced cDNA sequences with the entire genome sequence as a decoy.

#### Salmon parameters 

A benefit of using `salmon` is the ability to incorporate RNA-seq specific technical biases and correct counts accordingly. 
We chose to enable the [`--seqBias`](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias) and [`--gcBias`](https://salmon.readthedocs.io/en/latest/salmon.html#gcbias) flags, to correct for sequence-specific biases due to random hexamer primer and fragment-level GC biases, respectively. 
