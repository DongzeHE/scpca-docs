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

#### What genes are included in the reference transcriptome? 

The {ref}`reference transcriptome index <processing_information:reference transcriptome index>` that was used for alignment included 60,319 genes.
The reference transcriptome was constructed by extracting both spliced cDNA and intronic regions from the primary genome assembly ([see the code used to generate the reference transcriptome](https://github.com/AlexsLemonade/scpca-nf/blob/5dc49ffdcdb5c74dca86b0b79878f4d060029a53/bin/make_splici_fasta.R)).
In addition to protein-coding genes, this reference transcriptome index includes pseudogenes and non-coding RNA. 
Although there are 60,319 genes that are present in the reference transcriptome index, it is likely that many of these genes will not be detected in any given library.