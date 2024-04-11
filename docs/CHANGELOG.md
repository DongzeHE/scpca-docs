# CHANGELOG

As of November 2023, the `CHANGELOG` is a feature of our documentation we'll use to report and summarize changes to downloads from the ScPCA Portal.

You can find more information about how and when your download was prepared in the following places:

* The date your download was packaged (`Generated on: {date}`) is included at the top of the README in your download.
* The version of the [`AlexsLemonade/scpca-nf`](https://github.com/alexsLemonade/scpca-nf) pipeline used to process data in your download is included in the `workflow_version` column of the `single_cell_metadata.tsv` or `bulk_metadata.tsv` file in your download.
For more information about `AlexsLemonade/scpca-nf` versions, please see [the releases page on GitHub](https://github.com/AlexsLemonade/scpca-nf/releases).

<!-------------------------------------------------->
<!-- PUT THE NEW CHANGELOG ENTRY RIGHT BELOW THIS -->
<!-------------------------------------------------->


## 2024.04.11

* Cell type annotations are now included in each download.
Cells were annotated using both [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and [`CellAssign`](https://github.com/Irrationone/cellassign).
  * You can find more information about how cell types were annotated in the {ref}`cell type annotation procedures section on the Processing Information page<processing_information:cell type annotation>`.
  For more information on locating cell type annotations and any associated processing information in the downloaded objects see {ref}`the Single-cell gene expression file contents page<sce_file_conents:components of a singlecellexperiment object>`.
* Downloads also contain a separate cell type report providing more information about cell type annotations, including comparisons between different cell type annotations and diagnostic assessments of cell type annotation reliability.
* Sample metadata now includes two additional pieces of information which can be used to filter datasets: Whether the given sample is a patient-derived xenograft, and whether the sample is derived from a cell line.
* Processed `SingleCellExperiment` objects no longer include the full `miQC` result object in their metadata, but the `miQC` object is still available in the filtered `SingleCellExperiment` objects.
* Download files with `SingleCellExperiment` objects now use `bz2` compression.
This means the file sizes will be much smaller, but you may notice slower read times when loading them into R.
* This release additionally includes community-contributed projects.
Community-contributed projects are 10x Genomics single-cell or single-nuclei datasets that have been processed with the ScPCA pipeline.
Please refer to the [contributions page](https://scpca.alexslemonade.org/contribute) for more information about community contributions.


## 2024.03.08

* Downloads for most projects are now available in [`AnnData`](https://anndata.readthedocs.io/en/latest/index.html) format as HDF5 files.
Multiplexed samples are not yet supported.
* The sample metadata found in `single_cell_metadata.tsv` has been updated to include ontology term ids for age, sex, organism, ethnicity, diagnosis, and tissue location, when available.
See {ref}`the section describing Metadata on the Downloadable Files page<download_files:Metadata>`.
* All samples now have an assigned `participant_id`, which can be found in `single_cell_metadata.tsv`.
Previously, a `participant_id` was only assigned when multiple samples mapped to the same participant for most projects.
* All data files now include both the gene expression data and metadata for each sample (e.g., age, sex, organism, ethnicity, diagnosis, and tissue location).
For more information on the contents of the data files, see {ref}`the Single-cell gene expression file contents page<sce_file_contents:Single-cell gene expression file contents>`.
* Data files will include cell type annotations provided by submitters when applicable.

## 2023.11.10

* The README included in your download now contains the following:
	* More information about how to cite the ScPCA Portal (see also: {ref}`How to Cite <citation:how to cite>`).
	* The date your download was packaged (`Generated on: {date}`) at the top of the file.
* The root directory of your download will contain the date you accessed and downloaded data from the ScPCA Portal when uncompressed.
