# Downloadable files

The ScPCA Portal download packages include gene expression data, a QC report, and associated metadata for each processed sample.
Gene expression data is available as either [`SingleCellExperiment` objects (`.rds` files)](#singlecellexperiment-downloads) or [`AnnData` objects (`.hdf5` files)](#anndata-downloads).
These files are delivered as a zip file.
When you uncompress the zip file, the root directory name of your download will include the date you accessed the data on the ScPCA Portal.
We recommend you record this date in case there are future updates to the Portal that change the underlying data or if you need to cite the data in the future (see {ref}`How to Cite <citation:how to cite>` for more information).
Please see our {ref}`CHANGELOG <CHANGELOG:CHANGELOG>` for a summary of changes that impact downloads from the Portal.

For all downloads, sample folders (indicated by the `SCPCS` prefix) contain the files for all libraries (`SCPCL` prefix) derived from that biological sample.
Most samples only have one library that has been sequenced.
For [multiplexed sample libraries](#multiplexed-sample-libraries), the sample folder name will be an underscore-separated list of all samples found in the library files that the folder contains.
Note that multiplexed sample libraries are only available as `SingleCellExperiment` objects, and are not currently available as `AnnData` objects.

See the {ref}`FAQ section about samples and libraries <faq:What is the difference between samples and libraries?>` for more information.

The files shown below will be included with each library (example shown for a library with ID `SCPCL000000`):
- An unfiltered counts file: `SCPCL000000_unfiltered.rds` or `SCPCL00000_unfiltered_rna.hdf5`,
- A filtered counts file: `SCPCL000000_filtered.rds` or `SCPCL00000_filtered_rna.hdf5`,
- A processed counts file: `SCPCL000000_processed.rds` or `SCPCL00000_processed_rna.hdf5`,
- A quality control report: `SCPCL000000_qc.html`,
- A supplemental cell type report: `SCPCL000000_celltype-report.html`

Every download also includes a single `single_cell_metadata.tsv` file containing metadata for all libraries included in the download.

If downloading a project containing bulk RNA-seq data, two tab-separated value files, `bulk_quant.tsv` and `bulk_metadata.tsv`, will be included in the project download.
The `bulk_quant.tsv` file contains a gene by sample matrix (each row a gene, each column a sample) containing raw gene expression counts quantified by Salmon.
The `bulk_metadata.tsv` file contains associated metadata for all samples with bulk RNA-seq data.

See also {ref}`processing bulk RNA samples <processing_information:Bulk RNA samples>`.

The folder structure within the zip file is determined by whether individual samples or all samples associated with a project are selected for download.
Note that if a sample selected for download contains a spatial transcriptomics library, the files included will be different than pictured below.
See the [description of the Spatial transcriptomics output section below](#spatial-transcriptomics-libraries).

## `SingleCellExperiment` downloads

### Download folder structure for project downloads:
![project download folder](images/project-download-folder.png)

### Download folder structure for individual sample downloads:
![sample download folder](images/sample-download-folder.png)

## `AnnData` downloads

### Download folder structure for project downloads:
![project download folder](images/anndata-project-download-folder.png)

### Download folder structure for individual sample downloads:
![sample download folder](images/anndata-sample-download-folder.png)

### Download folder structure for individual sample downloads with CITE-seq (ADT) data:
![sample download folder](images/anndata-sample-citeseq-download-folder.png)

If downloading a sample that contains a CITE-seq library as an `AnnData` object (`hdf5` file), the quantified CITE-seq expression data is included as a separate file with the suffix `_adt.hdf5`.

## Gene expression data

Single-cell or single-nuclei gene expression data is provided as either [`SingleCellExperiment` objects (`.rds` files)](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) or [`AnnData` objects (`.hdf5` files)](https://anndata.readthedocs.io/en/latest/index.html).
Three files will be provided for each library included in the download - an unfiltered counts file, a filtered counts file, and a processed counts file.

The unfiltered counts file, `SCPCL000000_unfiltered.rds` or `SCPCL000000_unfiltered_rna.hdf5`, contains the counts matrix, where the rows correspond to genes or features and the columns correspond to cell barcodes.
Here, all potential cell barcodes that are identified after running {ref}`alevin-fry <processing_information:mapping and quantification using alevin-fry>` are included in the counts matrix.
The object also includes summary statistics for each cell barcode and gene, as well as metadata about that particular library, such as the reference index and software versions used for mapping and quantification.

The filtered counts file, `SCPCL000000_filtered.rds` or `SCPCL000000_filtered_rna.hdf5` contains a counts matrix with the same structure as above.
The cells in this file are those that remain after filtering using {ref}`emptyDrops <processing_information:filtering cells>`.
As a result, this file only contains cell barcodes that are likely to correspond to true cells.

The processed counts file, `SCPCL000000_processed.rds` or `SCPCL000000_processed_rna.hdf5`, contains both the raw and normalized counts matrices.
The filtered counts file is further filtered to remove low quality cells, such as those with a low number of genes detected or high mitochondrial content.
This file contains the raw and normalized counts data for cell barcodes that have passed both levels of filtering.
In addition to the counts matrices, the `SingleCellExperiment` or `AnnData` object stored in the file includes the results of dimensionality reduction using both principal component analysis (PCA) and UMAP.

See {ref}`Single-cell gene expression file contents <sce_file_contents:Single-cell gene expression file contents>` for more information about the contents of the `SingleCellExperiment` and `AnnData` objects and the included statistics and metadata.
See also {ref}`Using the provided RDS files in R <faq:how do i use the provided RDS files in r?>` and {ref}`Using the provided HDF5 files in Python <faq:how do i use the provided HDF5 files in python?>`.

## QC report

The included QC report, `SCPCL000000_qc.html`, serves as a general overview of each library, including processing information, summary statistics and general visualizations of cell metrics.

## Cell type report

The cell type report, `SCPCL000000_celltype-report.html`, includes an overview of cell type annotations present in the processed objects.
This report contains details on methodologies used for cell type annotation, information about reference sources, comparisons among cell type annotation methods, and diagnostic plots.
For more information on how cell types were annotated, see the section on {ref}`Cell type annotation <processing_information:cell type annotation>`.

## Metadata

The `single_cell_metadata.tsv` file is a tab-separated table with one row per library and the following columns.

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `scpca_sample_id` | Sample ID in the form `SCPCS000000`                            |
| `scpca_library_id` | Library ID in the form `SCPCL000000`                          |
| `seq_unit`        | `cell` for single-cell samples or `nucleus` for single-nuclei samples |
| `technology`      | 10x kit used to process library                                |
| `filtered_cell_count` | Number of cells after filtering with `emptyDrops`          |
| `submitter_id`    | Original sample identifier from submitter                      |
| `participant_id`  | Unique id corresponding to the donor from which the sample was obtained |
| `submitter`       | Submitter name/id                                              |
| `age`             | Age at time sample was obtained                                |
| `sex`             | Sex of patient that the sample was obtained from               |
| `diagnosis`       | Tumor type                                                     |
| `subdiagnosis`    | Subcategory of diagnosis or mutation status (if applicable)    |
| `tissue_location` | Where in the body the tumor sample was located                 |
| `disease_timing`  | At what stage of disease the sample was obtained, either diagnosis or recurrence |
| `organism`         | The organism the sample was obtained from (e.g., `Homo_sapiens`) |
| `development_stage_ontology_term_id` | [`HsapDv`](http://obofoundry.org/ontology/hsapdv.html) ontology term indicating the age at which the sample was collected. `unknown` indicates age is unavailable. |
| `sex_ontology_term_id`| [`PATO`](http://obofoundry.org/ontology/pato.html) term referring to the sex of the sample. `unknown` indicates sex is unavailable. |
| `organism_ontology_id`| NCBI taxonomy term for organism, e.g. [`NCBITaxon:9606`](https://ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_9606). |
| `self_reported_ethnicity_ontology_term_id` | For _Homo sapiens_ samples, a [`Hancestro` term](http://obofoundry.org/ontology/hancestro.html). `multiethnic` indicates more than one ethnicity is reported. `unknown` indicates unavailable ethnicity and `NA` is used for all other organisms. |
| `disease_ontology_term_id` | [`MONDO`](http://obofoundry.org/ontology/mondo.html) term indicating disease type. [`PATO:0000461`](https://ontobee.org/ontology/PATO?iri=http://purl.obolibrary.org/obo/PATO_0000461) is used for normal or healthy tissue. |
| `tissue_ontology_term_id` | [`UBERON`](http://obofoundry.org/ontology/uberon.html) term indicating tissue of origin. `NA` indicates tissue is unavailable.  |

Additional metadata may also be included, specific to the disease type and experimental design of the project.
Examples of this include treatment or outcome.
Metadata pertaining to processing will also be available in this table and inside of the `SingleCellExperiment` object.
See the {ref}`SingleCellExperiment experiment metadata <sce_file_contents:singlecellexperiment experiment metadata>` section for more information on metadata columns that can be found in the `SingleCellExperiment` object.
See the {ref}`AnnData experiment metadata <sce_file_contents:anndata experiment metadata>` section for more information on metadata columns that can be found in the `AnnData` object.

For projects with bulk RNA-seq data, the `bulk_metadata.tsv` file will be included for project downloads.
This file will contain fields equivalent to those found in the `single_cell_metadata.tsv` related to processing the sample, but will not contain patient or disease specific metadata (e.g. `age`, `sex`, `diagnosis`, `subdiagnosis`, `tissue_location`, or `disease_timing`).

## Multiplexed sample libraries

For libraries where multiple biological samples were combined via cellhashing or similar technology (see the {ref}`FAQ section about multiplexed samples <faq:What is a multiplexed sample?>`), the organization of the downloaded files and metadata is slightly different.
Note that multiplexed sample libraries are only available as `SingleCellExperiment` objects, and are not currently available as `AnnData` objects.

For project downloads, the counts and QC files will be organized by the _set_ of samples that comprise each library, rather than in individual sample folders.
These sample set folders are named with an underscore-separated list of the sample ids for the libraries within, _e.g._, `SCPCS999990_SCPCS999991_SCPCS999992`.
Bulk RNA-seq data, if present, will follow the [same format as bulk RNA-seq for single-sample libraries](#download-folder-structure-for-project-downloads).

![multiplexed project download folder](images/multiplexed-download-folder.png)

Because we do not perform demultiplexing to separate cells from multiplexed libraries into sample-specific count matrices, sample downloads from a project with multiplexed data will include all libraries that contain the sample of interest, but these libraries _will still contain cells from other samples_.

For more on the specific contents of multiplexed library `SingleCellExperiment` objects, see the {ref}`Additional SingleCellExperiment components for multiplexed libraries <sce_file_contents:additional singlecellexperiment components for multiplexed libraries>` section.

The [metadata file](#metadata) for multiplexed libraries (`single_cell_metadata.tsv`) will have the same format as for individual samples, but each row will represent a particular sample/library pair, meaning that there may be multiple rows for each `scpca_library_id`, one for each `scpca_sample_id` within that library.


## Merged object downloads

When downloading a full ScPCA project, you can choose to download each library as an individual file, or you can download {ref}`a single file containing all libraries merged into a single object<faq:INCOMING MERGED OBJECT SECTION>`. TODO!

Merged object downloads contain all single-cell or single-nuclei gene expression data for a given ScPCA project within a single object, provided as either a [`SingleCellExperiment` object (`.rds` file)](http://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html) or an [`AnnData` object (`.hdf5` file)](https://anndata.readthedocs.io/en/latest/index.html).

The object file, `SCPCP000000_merged.rds` or `SCPCP000000_merged_rna.hdf5`, contains both a raw and normalized counts matrix, each with combined counts for all libraries in an ScPCA project.
In addition to the counts matrices, the `SingleCellExperiment` or `AnnData` object stored in the file includes the results of library-weighted dimensionality reduction using both principal component analysis (PCA) and UMAP.
See the {ref}`section on merged object processing<processing_information:merged objects>` for more information about how merged objects were created.

If downloading a project that contains at least one CITE-seq library, the quantified CITE-seq expression data will also be merged.
In `SingleCellExperiment` objects (`rds` files), the  CITE-seq expression data is be provided as an alternative experiment in the same object as the gene expression data.
However, for `AnnData` objects, (`hdf5` files), the quantified CITE-seq expression is instead provided as a separate file called `SCPCP000000_merged_adt.hdf5`.

Every download includes a summary report, `SCPCL000000_merged-summary-report.html`, which provides a
 a brief summary of the libraries included in the merged object.
This includes a summary of the types of libraries (e.g., single-cell, single-nuclei, with CITE-seq) and sample diagnoses included in the object, as well as UMAP visualizations highlighting each library.

Every download also includes the individual [QC report](#qc-report) and, if applicable, [cell type annotation reports](#cell-type-report) for each library included in the merged object.

### Download folder structure for `SingleCellExperiment` merged downloads:
_image pending_

### Download folder structure for `AnnData` merged downloads:
_image pending_

### Download folder structure for `AnnData` merged downloads with CITE-seq (ADT) data:
_image pending_


### Merged object metadata

Similar to downloading the project with individual files for each sample, downloading the project as a merged object includes a single `single_cell_metadata.tsv` file containing metadata for all libraries included in the download.

If downloading a project containing bulk RNA-seq data, two tab-separated value files, `bulk_quant.tsv` and `bulk_metadata.tsv`, will be included in the project download.
The `bulk_quant.tsv` file contains a gene by sample matrix (each row a gene, each column a sample) containing raw gene expression counts quantified by Salmon.
The `bulk_metadata.tsv` file contains associated metadata for all samples with bulk RNA-seq data.


## Spatial transcriptomics libraries

If a sample includes a library processed using spatial transcriptomics, the spatial transcriptomics output files will be available as a separate download from the single-cell/single-nuclei gene expression data.

For all spatial transcriptomics libraries, a `SCPCL000000_spatial` folder will be nested inside the corresponding sample folder in the download.
Inside that folder will be the following folders and files:

- A `raw_feature_bc_matrix` folder containing the [unfiltered counts files](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices)
- A `filtered_feature_bc_matrix` folder containing the [filtered counts files](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices)
- A `spatial` folder containing [images and position information](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images)
- A `SCPCL000000_spaceranger_summary.html` file containing the [summary html report provided by Space Ranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/summary)
- A `SCPCL000000_metadata.json` file containing library processing information.

A full description of all files included in the download for spatial transcriptomics libraries can also be found in the [`spaceranger count` documentation](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count#outputs).

Every download also includes a single `spatial_metadata.tsv` file containing metadata for all libraries included in the download.

![sample download with spatial](images/spatial-download-folder.png)
