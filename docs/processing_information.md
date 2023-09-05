# Processing information

## Single-cell and single-nuclei RNA-seq

### Mapping and quantification using alevin-fry

We used [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/) to generate gene by cell counts matrices for all single-cell and single-nuclei samples.
In brief, we utilized [selective alignment](#selective-alignment) to the [`splici` index](#reference-transcriptome-index) for all single-cell and single-nuclei samples.

#### Reference transcriptome index

For all samples, we aligned FASTQ files to a reference transcriptome index referred to as the `splici` index.
The [`splici` index](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/) is built using transcripts from both spliced cDNA and intronic regions.
Inclusion of intronic regions in the index used for alignment allowed us to capture both reads from mature, spliced cDNA and nascent, unspliced cDNA.
Alignment of RNA-seq data to an index containing intronic regions has been shown to reduce spuriously detected genes ([He _et al._ 2021](https://doi.org/10.1101/2021.06.29.450377), [Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full#sec-5)).
In our hands, we have found that use of the `splici` index led to a more comparable distribution of unique genes found per cell to Cell Ranger than did use of an index obtained from spliced cDNA transcripts only.

#### Selective alignment

We mapped reads to the transcriptome index using `salmon` with the default "selective alignment" strategy.
Briefly, selective alignment uses a mapping score validated approach to identify maximal exact matches between reads and the provided index.
For all samples, we used selective alignment to the `splici` index.

A more detailed description of the mapping strategy invoked by `salmon` in conjunction with `alevin-fry` can be found in [Srivastava _et al._ (2020)](https://doi.org/10.1186/s13059-020-02151-8).

#### Alevin-fry parameters

After mapping FASTQ files using selective alignment to the `splici` index, we continued with the `alevin-fry` pipeline using the following parameters:

1. During the [`generate-permit-list` step of `alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/generate_permit_list.html), we used the `--unfiltered-pl` option, which returns any cell with at least 1 read found in a reference barcode list.
For our reference barcode list, we used a list of all possible cell barcodes from 10x Genomics.

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
To do this we used [`DropletUtils::emptyDropsCellRanger()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDropsCellRanger.html), a function that estimates the profile of cells containing ambient RNA and tests the likelihood of all other droplets as differing from the ambient profile ([Lun _et al._ 2019](https://doi.org/10.1186/s13059-019-1662-y)).
This function more closely mimics the filtering performed in Cell Ranger than does its predecessor [`DropletUtils::emptyDrops()`](https://www.bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets).
We consider droplets with an FDR less than or equal to 0.01 to be cell-containing droplets.
Only cells that pass this FDR threshold are included in the filtered counts matrix.

For some libraries, `DropletUtils::emptyDropsCellRanger()` may fail due to low numbers of droplets with reads or other violations of its assumptions.
For these libraries, only droplets containing at least 100 UMI are included in the filtered counts matrix.

### Processed gene expression data

In addition to the raw gene expression data, we also provide a processed `SingleCellExperiment` object with further filtering applied, a normalized counts matrix, and results from dimensionality reduction.

Prior to normalization, low-quality cells are removed from the gene by cell counts matrix.
To identify low-quality cells, we use [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html), a package that jointly models proportion of reads belonging to mitochondrial genes and number of unique genes detected.
Cells with a high likelihood of being compromised (greater than 0.75) and cells that do not pass a minimum number of unique genes detected threshold of 200 are removed from the counts matrix present in the processed `SingleCellExperiment` object.

Log-normalized counts are calculated using the deconvolution method presented in [Lun, Bach, and Marioni (2016)](https://doi.org/10.1186/s13059-016-0947-7).
The log-normalized counts are used to model variance of each gene prior to selecting the top 2000 highly variable genes (HVGs).
These HVGs are then used as input to principal component analysis, and the top 50 principal components are selected.
Finally, these principal components are used to calculate the [UMAP (Uniform Manifold Approximation and Projection)](http://bioconductor.org/books/3.13/OSCA.basic/dimensionality-reduction.html#uniform-manifold-approximation-and-projection) embeddings.

[Graph-based clustering](http://bioconductor.org/books/3.16/OSCA.basic/clustering.html#clustering-graph) is also performed, using the Louvain algorithm with 20 nearest neighbors and Jaccard weighting.

#### Cell type annotation

We perform cell type annotation with two complementary methods, where possible:

- [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html), a reference-based cell type annotation method
- [`CellAssign`](https://github.com/Irrationone/cellassign), a marker-gene--based cell type annotation method

For `SingleR` annotation, we identify an appropriate reference dataset from the [`celldex` package](http://bioconductor.org/packages/release/data/experiment/html/celldex.html) and train the classification model to use ontology IDs for annotation.
Cells which `SingleR` cannot confidently assign are labeled as `NA`.

For `CellAssign` annotation, we identify an appropriate set of marker genes for the given tissue type from the [`PanglaoDB`](https://panglaodb.se/) database.
We combine these marker genes with all "immune cell" `PanglaoDB` marker genes to create a full marker gene list for for cell type annotation.
During annotation, we additionally include an `"other"` cell type that does not express any of these marker genes.
As a consequence, cells which `CellAssign` cannot confidently annotate from the full marker gene list are labeled as `"other"`.
## ADT quantification from CITE-seq experiments

CITE-seq libraries with reads from antibody-derived tags (ADTs) were also quantified using  [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/), rounded to integer values.

Reference indices were constructed from the submitter-provided list of antibody barcode sequences corresponding to each library using the `--features` flag of `salmon index`.
Mapping to these indices followed the same procedures as for RNA-seq data, including mapping with [selective alignment](#selective-alignment) and subsequent [quantification via alevin-fry](#alevin-fry-parameters).

### Combining ADT counts with RNA counts

The unfiltered ADT and RNA-seq count matrices often include somewhat different sets of cell barcodes, due to stochastic variation in library construction and sequencing.
When normalizing these two count matrices to the same set of cells, we chose to prioritize RNA-seq results for broad comparability among libraries with and without ADT data.
Any cell barcodes that appeared only in ADT data were discarded.
Cell barcodes that were present only in the RNA-seq data (i.e., did not appear in the ADT data) were assigned zero counts for all ADTs.
When cells were [filtered based on RNA-seq content](#filtering-cells) after quantification, the ADT count matrix was filtered to match.

### Processed ADT data

An ambient profile representing antibody-derived tag (ADT) proportions present in the ambient solution is calculated from the unfiltered `SingleCellExperiment` object using [`DropletUtils::ambientProfileEmpty()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/ambientProfileEmpty.html).
Quality-control statistics were calculated with [`DropletUtils::cleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html) (with default parameters) using this ambient profile, along with negative/isotype control information, if present.
Low-quality cells identified by `DropletUtils::cleanTagCounts()` (those having high levels of ambient contamination or substantial negative/isotype control tags) are flagged but not removed except during normalization, as described below.
If `DropletUtils::cleanTagCounts()` cannot reliably determine which cells to filter, then no cells will be flagged for removal.

For all cells that would be retained if `DropletUtils::cleanTagCounts()` filtering were applied, log-normalized ADT counts are, by default, calculated using [median-based normalization](http://bioconductor.org/books/3.16/OSCA.advanced/integrating-with-protein-abundance.html#cite-seq-median-norm), again making use of the baseline ambient profile.
In order for this normalization to succeed, all median size factor values must be positive.
If any size factors are not positive or if ADT filtering failed, then only log-based normalization (with a pseudocount of one) will be performed.
Normalized counts for cells that would be filtered out by `DropletUtils::cleanTagCounts()` are assigned as `NA`.

## Multiplexed libraries

Multiplexed libraries, or libraries with cells or nuclei from more than one biological sample, were processed to allow demultiplexing using both hashtag oligonucleotide (HTO) results and genotype data.

### Hashtag oligonucleotide (HTO) quantification

HTO reads were also quantified using  [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) and [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/), rounded to integer values.
Reference indices were constructed from the submitter-provided list of HTO sequences corresponding to each library using the `--features` flag of `salmon index`.
Mapping to these indices followed the same procedures as for RNA-seq data, including mapping with [selective alignment](#selective-alignment) and subsequent [quantification via alevin-fry](#alevin-fry-parameters).

As with the [ADT data](#combining-adt-counts-with-rna-counts), we retained all cells with RNA-seq data, setting HTO counts to zero for any missing cell barcodes.
When cells were [filtered based on RNA-seq content](#filtering-cells) after quantification, the HTO count matrix was filtered to match.

### HTO demultiplexing

We performed HTO demultiplexing using both [`DropletUtils::hashedDrops()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html) and [`Seurat::HTODemux()`](https://rdrr.io/github/satijalab/seurat/man/HTODemux.html) only on the _filtered_ cells, using default parameters for each.

We report the demultiplexed sample calls and associated statistics for both algorithms, but do not separate the multiplexed library into individual samples.

### Genetic demultiplexing

For multiplex libraries where bulk RNA-seq data is available for the individual samples, we also performed demultiplexing analysis using genotype data following the methods described in [Weber _et al._ (2021)](https://doi.org/10.1093/gigascience/giab062):

- Bulk RNA-seq reads from each sample were mapped to the reference genome using `STAR` ([Dobin _et al._ 2012](https://doi.org/10.1093/bioinformatics/bts635))
- Variants among the samples within each pool were identified and genotyped with [`bcftools mpileup`](https://samtools.github.io/bcftools/bcftools.html#mpileup) ([Danecek _et al._ 2021](https://doi.org/10.1093/gigascience/giab008)) using the mapped bulk reads
- Pooled single-cell or single-nuclei RNA-seq reads were mapped to the reference genome using `STARsolo` ([Kaminow _et al._ 2021](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1))
- Individual cells were genotyped at the sites identified in the bulk RNA using [`cellsnp-lite`](https://cellsnp-lite.readthedocs.io) ([Huang _et al._ 2021](https://doi.org/10.1093/bioinformatics/btab358))
- Cell genotypes were used to call sample of origin with [`vireo`](https://vireosnp.readthedocs.io) ([Huang _et al._ 2019](https://doi.org/10.1186/s13059-019-1865-2))

The genetic demultiplexing calls are reported alongside HTO demultiplexing results for each library, but we again do not separate the individual samples.
For information on where the demultiplexing calls can be found, see {ref}`the section on demultiplexing results in the  SingleCellExperiment file contents. <sce_file_contents:demultiplexing results>`


## Spatial transcriptomics

### Mapping and quantification using Space Ranger

Processing spatial transcriptomics libraries requires two steps - gene expression quantification and tissue detection.
In the absence of independent tissue detection methods to use with Alevin-fry, we used [10x Genomics' Space Ranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger) to obtain both gene expression and spatial information.
[`spaceranger count`](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count) takes FASTQ files and a microscopic slide image as input and performs alignment, quantification, and tissue detection for each spot.
In contrast to `alevin-fry`, which maps reads to a [reference transcriptome index](#reference-transcriptome-index), Space Ranger aligns transcript reads to the reference genome using `STAR` ([Dobin _et al._ 2012](https://doi.org/10.1093/bioinformatics/bts635)).
See the 10x documentation for more information on how Space Ranger [quantifies gene expression](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/algorithms/overview) and [detects tissue](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/algorithms/imaging).

## Bulk RNA samples

### Preprocessing with fastp

Prior to quantifying gene expression for bulk RNA-seq samples, FASTQ files were pre-processed using [`fastp`](https://github.com/OpenGene/fastp) to perform adapter trimming, quality filtering, and length filtering.
For length filtering, trimmed reads shorter than 20 basepairs were removed by using the `--length_required 20` option.
All other filtering and trimming was performed using the default strategies enabled in `fastp`.

### Mapping and quantification using salmon

To quantify gene expression for bulk RNA-seq samples, we used [`salmon quant`](https://salmon.readthedocs.io/en/latest/salmon.html).
Here, we performed selective alignment of the trimmed and filtered FASTQ files to a decoy-aware reference transcriptome index ([Srivastava _et al._ 2020](https://doi.org/10.1186/s13059-020-02151-8)).
The [decoy-aware reference transcriptome](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode), was created from spliced cDNA sequences with the entire genome sequence as a decoy.

#### Salmon parameters

A benefit of using `salmon` is the ability to incorporate RNA-seq specific technical biases and correct counts accordingly.
We chose to enable the [`--seqBias`](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias) and [`--gcBias`](https://salmon.readthedocs.io/en/latest/salmon.html#gcbias) flags, to correct for sequence-specific biases due to random hexamer primer and fragment-level GC biases, respectively.
