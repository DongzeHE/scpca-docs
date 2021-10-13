# Frequently Asked Questions 

#### Why did we use Alevin-fry for processing? 

We aimed to process all of the data in the portal such that it is comparable to widely used pipelines, namely Cell Ranger from 10x Genomics.
In our own benchmarking, we found that [Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) produces very similar results to [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count), while allowing faster, more memory efficient processing of single-cell and single-nuclei RNA-sequencing data.
In the configuration that we are using ("selective alignment" mapping to a human transcriptome that includes introns), Alevin-fry uses approximately 12-16 GB of memory per sample and completes mapping and quantification in less than an hour. 
By contrast, Cell Ranger uses up to 25-30 GB of memory per sample and takes anywhere from 2-8 hours to align and quantify one sample.
Quantification of samples processed with both Alevin-fry and Cell Ranger resulted in similar distributions of mapped UMI count per cell and genes detected per cell for both tools.

![](https://github.com/AlexsLemonade/alsf-scpca/blob/dcebe62514b3a47a26b5923e7a1bdb1bb277e83b/analysis/docs-figures/plots/total_umi_per_cell.png?raw=true)

![](https://github.com/AlexsLemonade/alsf-scpca/blob/dcebe62514b3a47a26b5923e7a1bdb1bb277e83b/analysis/docs-figures/plots/total_genes_per_cell.png?raw=true)

We also compared the mean gene expression reported for each gene by both methods and observed a high correlation with a Pearson R correlation coefficient of 0.98.  

![](https://github.com/AlexsLemonade/alsf-scpca/blob/dcebe62514b3a47a26b5923e7a1bdb1bb277e83b/analysis/docs-figures/plots/gene_exp_correlation.png?raw=true)

Recent reports from others support our findings. 
[He _et al._ (2021)](https://doi.org/10.1101/2021.06.29.450377)) demonstrated that Alevin-fry can process single-cell and single-nuclei data more quickly and efficiently then other available methods, while also decreasing the false positive rate of gene detection that is commonly seen in methods that utilize transcriptome alignment.
