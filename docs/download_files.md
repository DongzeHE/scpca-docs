# Downloadable files

The ScPCA Portal download packages include gene expression data, a QC report, and associated metadata for each processed sample.
These files are delivered as a zip file.
The folder structure within the zip file is determined by whether individual samples or all samples associated with a project are selected for download.  

## Download folder structure for project downloads:
![docs-downloads-project](https://user-images.githubusercontent.com/15315514/143308420-a3cca10d-814f-4c52-b934-98d5e9cef1c5.png)

## Download folder structure for individual sample downloads:
![docs-downloads-sample](https://user-images.githubusercontent.com/15315514/143308420-a3cca10d-814f-4c52-b934-98d5e9cef1c5.png)

Sample folders (indicated by the `SCPCS` prefix) contain the files for all libraries (`SCPCL` prefix) derived from that biological sample.
Most samples only have one library that has been sequenced.
See the {ref}`FAQ section about samples and libraries <faq:What is the difference between samples and libraries?>` for more information.

The files associated with each library are (example shown for a library with ID `SCPCL000000`):
- An unfiltered counts file: `SCPCL000000_unfiltered.rds`,
- A filtered counts file: `SCPCL000000_filtered.rds`,
- A quality control report: `SCPCL000000_qc.html`,
- A metadata file: `SCPCL000000_metadata.json`.

Every download also includes a single `libraries_metadata.csv` file containing metadata for all libraries included in the download.

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

The `libraries_metadata.csv` file is a comma-separated table with one row per library and the following columns.

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

For each library, a `SCPCL000000_metadata.json` file that contains processing-associated metadata is also available.
Most fields can also be found in the metadata slot of the `SingleCellExperiment` objects stored in both the `SCPCL000000_unfiltered.rds` and `SCPCL000000_filtered.rds` files.
