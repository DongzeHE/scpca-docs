# CHANGELOG

As of November 2023, the `CHANGELOG` is a feature of our documentation we'll use to report and summarize changes to downloads from the ScPCA Portal.

You can find more information about how and when your download was prepared in the following places:

* The date your download was packaged (`Generated on: {date}`) is included at the top of the README in your download.
* The version of the [`AlexsLemonade/scpca-nf`](https://github.com/alexsLemonade/scpca-nf) pipeline used to process data in your download is included in the `workflow_version` column of the `single_cell_metadata.tsv` or `bulk_metadata.tsv` file in your download.
For more information about `AlexsLemonade/scpca-nf` versions, please see [the releases page on GitHub](https://github.com/AlexsLemonade/scpca-nf/releases).

<!-------------------------------------------------->
<!-- PUT THE NEW CHANGELOG ENTRY RIGHT BELOW THIS -->
<!-------------------------------------------------->

## PLACEHOLDER FOR RELEASE DATE

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
