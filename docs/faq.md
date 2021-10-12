# Frequently Asked Questions 

#### Why did we use Alevin-fry for processing? 

We have found that Alevin-fry offers a faster, more memory efficient approach to processing single-cell and single-nuclei RNA-sequencing data than using 10X Genomics' [Cellranger count,](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) without sacrificing accuracy. 
Alevin-fry uses approximately 12-16 GB of memory per sample and completes alignment and quantification in less than an hour, while `Cellranger count` uses up to 25-30 GB of memory per sample taking anywhere from 2-8 hours to align and quantify one sample.
Additionally, comparison of the quantification output from samples processed with both Alevin-fry and Cell Ranger revealed similar distributions of UMI/cell and genes detected/cell among both tools.

![](../images/total_umi_per_cell.png)

![](../images/total_genes_per_cell.png)

We also compared the mean gene expression reported for each gene by both methods and observed a high correlation with a Pearson R correlation coefficient of 0.98.  

*Insert Correlation graph comparing mean gene expression?*

Recent reports from others support our findings. 
[He _et al._ (2021)](https://doi.org/10.1101/2021.06.29.450377)) demonstrated that Alevin-fry can process single-cell and single-nuclei data more quickly and efficiently then other available methods, while also decreasing the false positive rate of gene detection that is commonly seen in methods that utilize transcriptome alignment.

