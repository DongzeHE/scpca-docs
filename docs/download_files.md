# Downloadable Files

The ScPCA portal download packages include gene expression data, a QC report, and associated metadata for each processed sample.
These files are delivered as a zip file.
The folder structure within the zip file is determined by whether a user selected to download all the samples associated with a project or selected to download individual samples.  

## Download folder structure for project downloads: 
_add in illustration of file structure for project_ 

## Download folder structure for individual sample downloads: 
_add in illustration of file structure for sample downloads_ 

Sample folders (indicated by the `SCPCS` prefix) contain the files for all libraries (`SCPCL` prefix) derived from that biological sample. 

The files associated with each library are:
- An unfiltered counts file: `SCPCL000000_unfiltered.rds`, 
- A filtered counts file: `SCPCL000000_filtered.rds`, 
- A quality control report: `SCPCL000000_qc.html`, 
- A metadata file: `SCPCL000000_metadata.json`. 
Every download also includes a `libraries_metadata.csv` file containing metadata associated with each library included in the download.

## Unfiltered Counts File 

The unfiltered counts matrix, `<scpca_library_id>_unfiltered.rds`, is delivered as a RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
Within the `SingleCellExperiment` is the counts matrix, where the rows correspond to genes or features and the columns correspond to cell barcodes. 
Here, all potential cell barcodes that are identified after running [Alevin-fry](processing_information.html/#alignment-and-quantification-using-alevin-fry) are included in the counts matrix. 
Summary statistics for each cell and gene can be found in the `colData` and `rowData` slots of the `SingleCellExperiment` object, respectively.
The `SingleCellExperiment` object also contains metadata about that particular library, including the versions of Salmon and Alevin-fry used for pre-processing, information about the index used for transcriptome alignment, and parameters used for Alevin-fry. 
See also [Using the provided RDS files in R.](link to FAQ)

## Filtered Counts Matrix

The filtered counts matrix is also delivered as a RDS file containing a `SingleCellExperiment` object.
Following filtering using [`emptyDrops`](processing_information.html/#filtering-cells), the filtered counts matrix and updated summary statistics for each cell and gene are output to the `<scpca_library_id>_filtered.rds` file.
As a result, this file only contains cell barcodes that are considered true cells.

## QC Report 

The included QC report serves as a general overview of each library, including processing information, summary statistics and general visualizations of cell metrics.

## Metadata

Metadata is provided in two ways. 
The `libraries_metadata.csv` file is a comma-separated table with one row per library and the following columns. 

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| scpca_sample_id | Sample ID in the form `SCPCS000000`                            |
| scpca_library_id | Library ID in the form `SCPCL000000`                          |
| submitter_id    | Original sample identifier from submitter                      |
| participant_id  | Original participant id, if there are multiple samples from the same participant                                                                          |
| submitter       | Submitter name/id                                              |
| age             | Age at time sample was obtained                                |
| sex             | Sex of patient sample is obtained from                      |
| diagnosis       | Tumor type                                                     |
| subdiagnosis    | Subcategory of diagnosis or mutation status (if applicable)   |
| tissue_location | Where tumor sample was located                                 |
| disease_timing  | What stage of disease was the sample obtained? At diagnosis or recurrence? |

Additional metadata may also be included, specific to the disease type and experimental design of the project. Examples of this include treatment or outcome. 

For each library, is also a `<scpca_library_id>_metadata.json` file that contains processing-associated metadata.
Most fields can also be found in the metadata slot of the `SingleCellExperiment` objects stored in both the `unfiltered.rds` and `filtered.rds` files.
