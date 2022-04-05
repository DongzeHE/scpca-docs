# Downloadable files

The ScPCA Portal download packages include gene expression data, a QC report, and associated metadata for each processed sample.
These files are delivered as a zip file.

For all downloads, sample folders (indicated by the `SCPCS` prefix) contain the files for all libraries (`SCPCL` prefix) derived from that biological sample.
Most samples only have one library that has been sequenced.
See the {ref}`FAQ section about samples and libraries <faq:What is the difference between samples and libraries?>` for more information.

The files associated with each library are (example shown for a library with ID `SCPCL000000`):
- An unfiltered counts file: `SCPCL000000_unfiltered.rds`,
- A filtered counts file: `SCPCL000000_filtered.rds`,
- A quality control report: `SCPCL000000_qc.html`,

Every download also includes a single `single_cell_metadata.tsv` file containing metadata for all libraries included in the download.

If a sample includes a library that was processed using spatial transcriptomics, there will be an additional folder in the sample directory, `SCPCL000000_spatial`. 
For more information about the contents of this folder, see the [spatial transcriptomics libraries section below](#spatial-transcriptomics-libraries). 

The folder structure within the zip file is determined by whether individual samples or all samples associated with a project are selected for download.  

## Download folder structure for project downloads:
![docs-downloads-project](https://user-images.githubusercontent.com/15315514/156806533-f45e1bef-9a96-490d-b9ad-5aabec2d28a2.png)

If a project contains bulk RNA-seq data, two tab-separated value files, `bulk_quant.tsv` and `bulk_metadata.tsv`, will be included in the download. 
The `bulk_quant.tsv` file contains a gene by sample matrix (each row a gene, each column a sample) containing raw gene expression counts quantified by Salmon.
The `bulk_metadata.tsv` file contains associated metadata for all samples with bulk RNA-seq data.  
 
See also {ref}`processing bulk RNA samples <processing_information:Bulk RNA samples>`.   

## Download folder structure for individual sample downloads:
![docs-downloads-sample](https://user-images.githubusercontent.com/15315514/156806534-7ef8d56b-f8f9-4458-b3f5-3f0f42cb7ede.png)

Note that if a sample selected for download contains a spatial transcriptomics library, the files included will be different than pictured above. 
See the [description of the Spatial transcriptomics output section below](#download-contents-for-samples-with-spatial-transcriptomics). 

## Gene expression data

Single-cell or single-nuclei gene expression data is provided in two forms - as an unfiltered counts file and a filtered counts file.

The unfiltered counts file, `SCPCL000000_unfiltered.rds`, is an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
Within the `SingleCellExperiment` object is the counts matrix, where the rows correspond to genes or features and the columns correspond to cell barcodes.
Here, all potential cell barcodes that are identified after running {ref}`Alevin-fry <processing_information:mapping and quantification using alevin-fry>` are included in the counts matrix.
The object also includes summary statistics for each cell barcode and gene, as well as metadata about that particular library, such as the reference index and software versions used for mapping and quantification.

The filtered counts file, `SCPCL000000_filtered.rds` is also an RDS file containing a `SingleCellExperiment` object with the same structure as above.
The cells in this file are those that remain after filtering using {ref}`emptyDrops <processing_information:filtering cells>`.
As a result, this file only contains cell barcodes that are likely to correspond to true cells.

See {ref}`Single-cell gene expression file contents <sce_file_contents:Single-cell gene expression file contents>` for more information about the contents of the `SingleCellExperiment` objects and the included statistics and metadata.
See also {ref}`Using the provided RDS files in R <faq:how do i use the provided RDS files in r?>`.

## QC Report

The included QC report serves as a general overview of each library, including processing information, summary statistics and general visualizations of cell metrics.

## Metadata

The `single_cell_metadata.tsv` file is a comma-separated table with one row per library and the following columns.

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| scpca_sample_id | Sample ID in the form `SCPCS000000`                            |
| scpca_library_id | Library ID in the form `SCPCL000000`                          |
| seq_unit        | `cell` for single-cell samples or `nucleus` for single-nucleus samples |
| technology      | 10X kit used to process library                                |
| filtered_cell_count | Number of cells after filtering with `emptyDrops`          |
| submitter_id    | Original sample identifier from submitter                      |
| participant_id  | Original participant id, if there are multiple samples from the same participant                                                                        |
| submitter       | Submitter name/id                                              |
| age             | Age at time sample was obtained                                |
| sex             | Sex of patient that the sample was obtained from               |
| diagnosis       | Tumor type                                                     |
| subdiagnosis    | Subcategory of diagnosis or mutation status (if applicable)    |
| tissue_location | Where in the body the tumor sample was located                 |
| disease_timing  | What stage of disease was the sample obtained? At diagnosis or recurrence? |

Additional metadata may also be included, specific to the disease type and experimental design of the project.
Examples of this include treatment or outcome.
Metadata pertaining to processing will also be available in this table and inside of the `SingleCellExperiment` object.
See the {ref}`Experiment metadata <sce_file_contents:experiment metadata>` section for more information on metadata columns that can be found in this file as well as inside the `SingleCellExperiment` object.

For projects with bulk RNA-seq data, the `bulk_metadata.tsv` file will be included for project downloads. 
This file will contain fields equivalent to those found in the `single_cell_metadata.tsv` related to processing the sample, but will not contain patient or disease specific metadata (e.g. `age`, `sex`, `diagnosis`, `subdiagnosis`, `tissue_location`, or `disease_timing`).

## Spatial transcriptomics libraries

For all spatial transcriptomics libraries, a `SCPCL000000_spatial` folder will be nested inside the corresponding sample folder in the download. 
Inside that folder will be the following folders and files: 

- A `raw_feature_bc_matrix` folder containing the [unfiltered counts files](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices)
- A `filtered_feature_bc_matrix` folder containing the [filtered counts files](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices)
- A `spatial` folder containing [images and position information](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images)
- A `SCPCL000000_spaceranger_summary.html` file containing the Space Ranger [summary html report](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/summary)

A full description of each of files included in the download for spatial transcriptomics libraries can also be found in the [`spaceranger count` documentation](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count#outputs). 

_Coming Soon: Illustration of example download with spatial library_