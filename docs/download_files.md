# Downloadable Files

Users will be able to download the gene expression data, QC report, and associated metadata for all samples that have been processed.
These files are delivered as a zip file.
The folder structure within the zip file is determined by whether a user selected to download all the samples associated with a project or selected to download individual samples.  

## Download folder structure for project downloads: 
_add in illustration of file structure for project_ 

## Download folder structure for individual sample downloads: 
_add in illustration of file structure for sample downloads_ 

Folders for each sample (indicated by the `SCPCS` prefix) will contain the files for all libraries corresponding to that sample. See also [difference between samples and libraries.](link to FAQ)
Each library will include the following files: `<scpca_library_id>_unfiltered.rds`, `<scpca_library_id>_filtered.rds`, `<scpca_library_id>_qc.html`, and `<scpca_library_id>_metadata.json`. 
Every download will also include a `libraries_metadata.csv` file containing metadata associated with each library included in the download.

### Unfiltered Counts Matrix 

The unfiltered counts matrix is delivered as a RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html).
Within the `SingleCellExperiment` is the counts matrix, where the rows correspond to genes or features and the columns correspond to cell barcodes. 
Here, all potential cell barcodes that are identifed after running [Alevin-fry](processing_information.md/#alignment-and-quantification-using-alevin-fry) are included in the counts matrix. 
Summary statistics for each cell and gene can be found in the `colData` and `rowData` of the `SingleCellExperiment` object.
The `SingleCellExperiment` object also contains metadata about that particular library, including the versions of Salmon and Alevin-fry used for pre-processing, information about the index used for transcriptome alignment, and parameters used for Alevin-fry. 
See also [Using the provided RDS files in R.](link to FAQ)

### Filtered Counts Matrix

The filtered counts matrix is also delivered as a rds file containing a `SingleCellExperiment` object.
Following filtering using [`emptyDrops`](processing_information.md/#filtering-cells), the filtered counts matrix and updated summary statistics for each cell and gene are output to the `<scpca_library_id>_filtered.rds` file.
As a result, this file only contains cell barcodes that are considered true cells.

### QC Report 

### Metadata
