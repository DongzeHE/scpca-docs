# Merged objects

Each merged object contains _combined information_ for all individual samples in the given ScPCA project.
While each individual object, as described on the {ref}`Single-cell gene expression file contents page <sce_file_contents:Single-cell gene expression file contents>`, contains quantified gene expression results for a single library, each merged object contains all gene expression results, including gene expression counts and metadata, for all libraries and samples in the given ScPCA project.
This information includes quantified gene expression data, cell and gene metrics, and associated metadata for all libraries.
See {ref}`the section on merged object processing <processing_information:merged objects>` for more information on how these objects were prepared.


Merged objects are provided in two formats:
  - As an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.17/OSCA.intro/the-singlecellexperiment-class.html) for use in R.
  - As an H5AD file containing an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html) for use in Python.


Below we present some details about the specific contents of the objects we provide.

## Components of a SingleCellExperiment merged object

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
merged_sce <- readRDS("SCPCP000000_merged.rds")
```

### SingleCellExperiment expression counts

Merged `SingleCellExperiment` objects contain two main assays, `counts` and `logcounts`, each containing RNA-seq expression data for all libraries in the given ScPCA project combined into a single matrix.
The `counts` assay contains the primary raw counts represented as integers, and the `logcounts` assay contains normalized counts as described in {ref}`the data post-processing section <processing_information:processed gene expression data>`.

The `counts` assay includes reads aligned to both spliced and unspliced cDNA (see the section on {ref}`Post Alevin-fry processing <processing_information:post alevin-fry processing>`).
Each assay is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `counts` and `logcounts` assays can be accessed with the following R code:

```r
counts(merged_sce) # combined counts matrix
logcounts(merged_sce) # combined logcounts matrix
```

Column names are cell barcode sequences prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed with the following R code:

```r
colnames(merged_sce) # matrix column names
rownames(merged_sce) # matrix row names
```

There is also a `spliced` assay which contains the counts matrix with only reads from spliced cDNA.

### SingleCellExperiment cell metrics


Cell metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `colData` slot, where row names are the cell barcode prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`.
This `DataFrame` also contains additional sample metadata information stored in the `colData` slot for all projects that do not contain multiplexed libraries.
Read more about the included sample metadata in the [`Sample metadata` section](#singlecellexperiment-sample-metadata),

```r
colData(merged_sce) # cell metrics
```

The following per-cell data columns are included for each cell.
Columns representing quality control statistics were calculated using the [`scuttle::addPerCellQCMetrics()`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) function.

| Column name                      | Contents                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `barcodes`                       | The cell barcode                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `library_id`                     | Library ID in the form `SCPCL000000`                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `cell_id`                        | Unique ID for each cell in the format `SCPCL000000-{barcode}`. This value matches the `colData` row names                                                                                                                                                                                                                                                                                                                                              |
| `sum`                            | UMI count for RNA-seq data                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `detected`                       | Number of genes detected (gene count > 0 )                                                                                                                                                                                                                                                                                                                                                                                                             |
| `subsets_mito_sum`               | UMI count of mitochondrial genes                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `subsets_mito_detected`          | Number of mitochondrial genes detected                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `subsets_mito_percent`           | Percent of all UMI counts assigned to mitochondrial genes                                                                                                                                                                                                                                                                                                                                                                                              |
| `total`                          | Total UMI count for RNA-seq data and any alternative experiments (i.e., ADT data from CITE-seq)                                                                                                                                                                                                                                                                                                                                                        |
| `prob_compromised`               | Probability that a cell is compromised (i.e., dead or damaged), as calculated by `miQC`                                                                                                                                                                                                                                                                                                                                                                |
| `miQC_pass`                      | Indicates whether the cell passed the default miQC filtering. `TRUE` is assigned to cells with a low probability of being compromised (`prob_compromised` < 0.75) or [sufficiently low mitochondrial content](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html#preventing-exclusion-of-low-mito-cells)                                                                                                                 |
| `scpca_filter`                   | Labels cells as either `Keep` or `Remove` based on filtering criteria (`prob_compromised` < 0.75 and number of unique genes detected > 200)                                                                                                                                                                                                                                                                                                            |
| `adt_scpca_filter`               | If CITE-seq was performed, labels cells as either `Keep` or `Remove` based on ADT filtering criteria (`discard = TRUE` as determined by [`DropletUtils::CleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html))                                                                                                                                                                                                    |

Unlike for {ref}`individual SCE objects<sce_file_contents:singlecellexperiment cell metrics>`, cluster assignments are not included in the `colData`.


### SingleCellExperiment gene information and metrics

Gene information and metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `rowData` slot, with the Ensembl ID as the names of the rows.

```r
rowData(merged_sce) # gene metrics
```

The following columns are included for all genes.
The columns `mean` and `detected` will appear for each library ID included in the merged object, named as shown in the table below.
However, there will only be a single `gene_symbol` and `gene_ids` column, as this information equally pertains to all libraries.
Metrics were calculated for each library using the [`scuttle::addPerFeatureQCMetrics`](https://rdrr.io/github/LTLA/scuttle/man/addPerFeatureQCMetrics.html) function.

| Column name            | Contents                                                                                    |
| ---------------------- | ------------------------------------------------------------------------------------------- |
| `gene_ids`             | Ensembl gene ID                                                                             |
| `gene_symbol`          | [HUGO](https://www.genenames.org) gene symbol, if defined                                   |
| `SCPCL000000-mean`     | Mean count across all cells/droplets for library `SCPCL000000`                              |
| `SCPCL000000-detected` | Percent of cells in which the gene was detected (gene count > 0 ) for library `SCPCL000000` |


### SingleCellExperiment experiment metadata

```r
metadata(merged_sce) # experiment metadata
```

| Item name                      | Contents                                                                                                                                                                                                                                                                                                                                               |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `library_id`                   | A vector of library IDs which are included in the merged object, each in the form `SCPCL000000`                                                                                                                                                                                                                                                        |
| `sample_id`                    | A vector of sample IDs which are included in the merged object, each in the form  `SCPCS000000`                                                                                                                                                                                                                                                       |
| `library_metadata`             | A list of the library metadata for each library. Each list is named with the appropriate library ID and contains the library metadata fields for the given library {ref}`as they would appear in an individual library object<sce_file_contents:singlecellexperiment experiment metadata>`. See the table below for a full description of its contents |
| `sample_metadata`              | A data frame of additional sample metadata, as described in the [`Sample metadata` section below](#singlecellexperiment-sample-metadata). Only present for projects with multiplexed libraries                                                                                                                                                         |
| `merged_highly_variable_genes` | A vector of highly variable genes used for performing dimensionality reduction on the merged object, determined using `scran::modelGeneVar`, specifying each library as a separate block, and `scran::getTopHVGs`                                                                                                                                      |

To access the `library_metadata` field for a specific library, use the following code:

```r
# Access individual library metadata for SCPCL000000
metadata(merged_sce)$library_metadata$SCPCL000000
```

Each such list will contain the following fields:


| Item name                      | Contents                                                                                                                                                                                                                                                                                                                                                                           |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample_id`                    | Sample ID in the form `SCPCS000000`                                                                                                                                                                                                                                                                                                                                                |
| `library_id`                   | Library ID in the form `SCPCL000000`                                                                                                                                                                                                                                                                                                                                               |
| `project_id`                   | Project ID in the form `SCPCP000000`                                                                                                                                                                                                                                                                                                                                               |
| `salmon_version`               | Version of `salmon` used for initial mapping                                                                                                                                                                                                                                                                                                                                       |
| `reference_index`              | Transcriptome reference file used for mapping                                                                                                                                                                                                                                                                                                                                      |
| `total_reads`                  | Total number of reads processed by `salmon`                                                                                                                                                                                                                                                                                                                                        |
| `mapped_reads`                 | Number of reads successfully mapped                                                                                                                                                                                                                                                                                                                                                |
| `mapping_tool`                 | Pipeline used for mapping and quantification (`alevin-fry` for all current data in ScPCA)                                                                                                                                                                                                                                                                                          |
| `alevinfry_version`            | Version of `alevin-fry` used for mapping and quantification                                                                                                                                                                                                                                                                                                                        |
| `af_permit_type`               | `alevin-fry generate-permit-list` method used for filtering cell barcodes                                                                                                                                                                                                                                                                                                          |
| `af_resolution`                | `alevin-fry quant` resolution mode used                                                                                                                                                                                                                                                                                                                                            |
| `usa_mode`                     | Boolean indicating whether quantification was done using `alevin-fry` USA mode                                                                                                                                                                                                                                                                                                     |
| `af_num_cells`                 | Number of cells reported by `alevin-fry`                                                                                                                                                                                                                                                                                                                                           |
| `tech_version`                 | A string indicating the technology and version used for the single-cell library, such as 10Xv2, 10Xv3, or 10Xv3.1                                                                                                                                                                                                                                                                  |
| `assay_ontology_term_id`       | A string indicating the [Experimental Factor Ontology](https://www.ebi.ac.uk/ols/ontologies/efo) term ID associated with the `tech_version`                                                                                                                                                                                                                                        |
| `seq_unit`                     | `cell` for single-cell samples or `nucleus` for single-nucleus samples                                                                                                                                                                                                                                                                                                             |
| `transcript_type`              | Transcripts included in gene counts: `spliced` for single-cell samples and `unspliced` for single-nuclei                                                                                                                                                                                                                                                                           |
| `sample_type`                  | A string indicating the type of sample, with one of the following values: `"patient-derived xenograft"`, `"cell line"`, or `"patient tissue"`. If the library is multiplexed, this will be a named vector giving the sample type for each sample ID in the library. A value of `"Not provided"` indicates that this information is not available                                                                                                                 |
| `filtering_method`             | The method used for cell filtering. One of `emptyDrops`, `emptyDropsCellRanger`, or `UMI cutoff`                                                                                                                                                                                                                                                                                   |
| `umi_cutoff`                   | The minimum UMI count per cell used as a threshold for removing empty droplets. Only present for objects where the `filtering_method` is `UMI cutoff`                                                                                                                                                                                                                              |
| `prob_compromised_cutoff`      | The minimum cutoff for the probability of a cell being compromised, as calculated by `miQC`                                                                                                                                                                                                                                                                                        |
| `scpca_filter_method`          | Method used by the Data Lab to filter low quality cells prior to normalization. Either `miQC` or `Minimum_gene_cutoff`                                                                                                                                                                                                                                                             |
| `adt_scpca_filter_method`      | If CITE-seq was performed, the method used by the Data Lab to identify cells to be filtered prior to normalization, based on ADT counts. Either `cleanTagCounts with isotype controls` or `cleanTagCounts without isotype controls`. If filtering failed (i.e. `DropletUtils::cleanTagCounts()` could not reliably determine which cells to filter), the value will be `No filter` |
| `min_gene_cutoff`              | The minimum cutoff for the number of unique genes detected per cell                                                                                                                                                                                                                                                                                                                |
| `normalization`                | The method used for normalization of raw RNA counts. Either `deconvolution`, described in [Lun, Bach, and Marioni (2016)](https://doi.org/10.1186/s13059-016-0947-7), or `log-normalization`                                                                                                                                                                                       |
| `adt_normalization`            | If CITE-seq was performed, the method used for normalization of raw ADT counts. Either `median-based` or `log-normalization`, as explained in the {ref}`processed ADT data section <processing_information:Processed ADT data>`                                                                                                                                                    |
| `highly_variable_genes`        | A list of highly variable genes used for dimensionality reduction, determined using `scran::modelGeneVar` and `scran::getTopHVGs`                                                                                                                                                                                                                                                  |
| `celltype_methods`             | If cell type annotation was performed, a vector of the methods used for annotation. May include `"submitter"`, `"singler"` and/or `"cellassign"`                                                                                                                                                                                                                                   |
| `singler_results`              | If cell typing with `SingleR` was performed, the full result object returned by `SingleR` annotation                                                                                                                                                                                                                                                                               |
| `singler_reference`            | If cell typing with `SingleR` was performed, the name of the reference dataset used for annotation                                                                                                                                                                                                                                                                                 |
| `singler_reference_label`      | If cell typing with `SingleR` was performed, the name of the label in the reference dataset used for annotation                                                                                                                                                                                                                                                                    |
| `singler_reference_source`     | If cell typing with `SingleR` was performed, the source of the reference dataset (default is [`celldex`](http://bioconductor.org/packages/release/data/experiment/html/celldex.html))                                                                                                                                                                                              |
| `singler_reference_version`    | If cell typing with `SingleR` was performed, the version of `celldex` used to create the reference dataset source, with periods replaced as dashes (`-`)                                                                                                                                                                                                                           |
| `cellassign_predictions`       | If cell typing with `CellAssign` was performed and completed successfully, the full matrix of predictions across cells and cell types                                                                                                                                                                                                                                              |
| `cellassign_reference`         | If cell typing with `CellAssign` was performed and completed successfully, the reference name as established by the Data Lab used for cell type annotation                                                                                                                                                                                                                         |
| `cellassign_reference_organs`  | If cell typing with `CellAssign` was performed and completed successfully, a comma-separated list of organs and/or tissue compartments from which marker genes were obtained to create the reference                                                                                                                                                                               |
| `cellassign_reference_source`  | If cell typing with `CellAssign` was performed and completed successfully, the source of the reference dataset (default is [`PanglaoDB`](https://panglaodb.se/))                                                                                                                                                                                                                   |
| `cellassign_reference_version` | If cell typing with `CellAssign` was performed and completed successfully, the version of the reference dataset source. For references obtained from `PanglaoDB`, the version scheme is a date in ISO8601 format                                                                                                                                                                   |

Unlike for {ref}`individual SingleCellExperiment objects<sce_file_contents:singlecellexperiment sample metadata>`, cluster algorithm parameters are not included in these metadata lists because clusters themselves are not included in the merged object.


### SingleCellExperiment sample metadata

Sample metadata describing each sample included in the merged object is stored in one of two locations, depending on the sample type:

If the project does_not_ contain multiplexed libraries, this information can be found as additional columns in the `colData` slot's `DataFrame`, along with [cell metrics](#singlecellexperiment-cell-metrics).

```r
colData(merged_sce) # sample metadata only for projects without multiplexing
```

If the project contains multiplexed libraries, this information is stored in the `metadata` slot in the `sample_metadata` field as a `data.frame`.
Similar to merged objects with multiplexed libraries, all {ref} `individual library objects<sce_file_contents:singlecellexperiment sample metadata>` will contain this sample metadata in the `SingleCellExperiment` object's `metadata` slot.

```r
metadata(merged_sce)$sample_metadata # sample metadata only for projects with multiplexed samples
```


| Column name                                | Contents                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample_id`                                | Sample ID in the form `SCPCS000000`                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `scpca_project_id`                         | Project ID in the form `SCPCP000000`                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `tech_version`                             | A string indicating the technology and version used for the sample's single-cell library, such as 10Xv2, 10Xv3, or 10Xv3.1                                                                                                                                                                                                                                                                                                                             |
| `assay_ontology_term_id`                   | A string indicating the [Experimental Factor Ontology](https://www.ebi.ac.uk/ols/ontologies/efo) term ID associated with the `tech_version`                                                                                                                                                                                                                                                                                                            |
| `suspension_type`                          | `cell` for single-cell samples or `nucleus` for single-nucleus samples                                                                                                                                                                                                                                                                                                                                                                                 |
| `additional_modalities`                    | Any additional modalities associated with the library, represented as alternative experiment names such as `"adt"` or `"cellhash"`. If there are no additional modalities, this value will be `NA`. If the library has multiple additional modalities, this will be a single string with modalities separated by a semicolon, e.g. `"adt;cellhash"`                                                                                                    |
| `participant_id`                           | Unique ID corresponding to the donor from which the sample was obtained                                                                                                                                                                                                                                                                                                                                                                                |
| `submitter_id`                             | Original sample identifier from submitter                                                                                                                                                                                                                                                                                                                                                                                                              |
| `submitter`                                | Submitter name/ID                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `age`                                      | Age at time sample was obtained                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `sex`                                      | Sex of patient that the sample was obtained from                                                                                                                                                                                                                                                                                                                                                                                                       |
| `diagnosis`                                | Tumor type                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `subdiagnosis`                             | Subcategory of diagnosis or mutation status (if applicable)                                                                                                                                                                                                                                                                                                                                                                                            |
| `tissue_location`                          | Where in the body the tumor sample was located                                                                                                                                                                                                                                                                                                                                                                                                         |
| `disease_timing`                           | At what stage of disease the sample was obtained, either diagnosis or recurrence                                                                                                                                                                                                                                                                                                                                                                       |
| `organism`                                 | The organism the sample was obtained from (e.g., `Homo_sapiens`)                                                                                                                                                                                                                                                                                                                                                                                       |
| `is_xenograft`                             | Whether the sample is a patient-derived xenograft                                                                                                                                                                                                                                                                                                                                                                                                      |
| `is_cell_line`                             | Whether the sample was derived from a cell line                                                                                                                                                                                                                                                                                                                                                                                                        |
| `development_stage_ontology_term_id`       | [`HsapDv` ontology](http://obofoundry.org/ontology/hsapdv.html) term indicating developmental stage. If unavailable, `unknown` is used                                                                                                                                                                                                                                                                                                                 |
| `sex_ontology_term_id`                     | [`PATO`](http://obofoundry.org/ontology/pato.html) term referring to the sex of the sample. If unavailable, `unknown` is used                                                                                                                                                                                                                                                                                                                          |
| `organism_ontology_id`                     | [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) term for organism, e.g. [`NCBITaxon:9606`](http://purl.obolibrary.org/obo/NCBITaxon_9606)                                                                                                                                                                                                                                                                                                       |
| `self_reported_ethnicity_ontology_term_id` | For _Homo sapiens_, a [`Hancestro` term](http://obofoundry.org/ontology/hancestro.html). `multiethnic` indicates more than one ethnicity is reported. `unknown` indicates unavailable ethnicity and `NA` is used for all other organisms                                                                                                                                                                                                               |
| `disease_ontology_term_id`                 | [`MONDO`](http://obofoundry.org/ontology/mondo.html) term indicating disease type. [`PATO:0000461`](http://purl.obolibrary.org/obo/PATO_0000461) indicates normal or healthy tissue. If unavailable, `NA` is used                                                                                                                                                                                                                                      |
| `tissue_ontology_term_id`                  | [`UBERON`](http://obofoundry.org/ontology/uberon.html) term indicating tissue of origin. If unavailable, `NA` is used                                                                                                                                                                                                                                                                                                                                  |
| `submitter_celltype_annotation`            | If available, cell type annotations obtained from the group that submitted the original data. Cells that the submitter did not annotate are labeled as `"Submitter-excluded"`                                                                                                                                                                                                                                                                          |
| `singler_celltype_annotation`              | If cell typing with `SingleR` was performed, the annotated cell type. Cells labeled as `NA` are those which `SingleR` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"`                                                                                                                           |
| `singler_celltype_ontology`                | If cell typing with `SingleR` was performed with ontology labels, the annotated cell type's ontology ID. Cells labeled as `NA` are those which `SingleR` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"`                                                                                        |
| `cellassign_celltype_annotation`           | If cell typing with `CellAssign` was performed, the annotated cell type. Cells labeled as `"other"` are those which `CellAssign` could not confidently annotate. If cell typing was performed for some but not all libraries in the merged object, libraries without annotations will be labeled `"Cell type annotation not performed"`. If `CellAssign` was unable to complete successfully for a given library, cells will be labeled as `"Not run"` |
| `cellassign_max_prediction`                | If cell typing with `CellAssign` was performed and completed successfully, the annotation's prediction score (probability)                                                                                                                                                                                                                                                                                                                             |

Some merged objects may have some additional sample metadata columns specific to the given ScPCA project's disease type and experimental design.
Examples of this include treatment or outcome.


### SingleCellExperiment dimensionality reduction results

The `reducedDim` slot of the merged object will contain both principal component analysis (`PCA`) and `UMAP` results.

PCA results were calculated using [`batchelor::multiBatchPCA`](https://rdrr.io/bioc/batchelor/man/multiBatchPCA.html), specifying libraries as batches to ensure that each library in the merged object was equally weighted, and specifying a list of highly variable genes.
The highly variable genes were selected in a library-aware manner with `scran::modelGeneVar` and `scran::getTopHVGs`.
The vector of highly variable genes are stored in the `SingleCellExperiment` object in `metadata(merged_sce)$merged_highly_variable_genes`.
The following command can be used to access the PCA results:

```r
reducedDim(merged_sce, "PCA")
```

UMAP results were calculated using `scater::runUMAP()`, with the batch-aware PCA results as input rather than the full gene expression matrix.
The following command can be used to access the UMAP results:

```r
reducedDim(merged_sce, "UMAP")
```

### Additional SingleCellExperiment components for CITE-seq libraries (with ADT tags)

ADT data from CITE-seq experiments, when present, is included within the merged `SingleCellExperiment` object as an "Alternative Experiment" named `"adt"`, which can be accessed with the following command:

```r
altExp(merged_sce, "adt") # adt experiment
```

Within this, the primary expression matrix is again found in the `counts` assay, and the normalized expression matrix is found in the `logcounts` assay.
Each expression matrix contains the ADT expression data from the CITE-seq experiment for all libraries in the given ScPCA project combined into a single matrix.
For each assay, each column corresponds to a cell or droplet (in the same order as the parent `SingleCellExperiment`) and each row corresponds to an antibody derived tag (ADT).
Column names are again cell barcode sequences prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`, and row names are the antibody targets for each ADT.
In some cases, only some libraries in the merged object will have associated CITE-seq data.
Libraries which do not have CITE-seq expression data are included in the assay matrices with all `NA` count values.

Only cells which are denoted as "Keep" in the `colData(merged_sce)$adt_scpca_filter` column (as described [above](#singlecellexperiment-cell-metrics)) have normalized expression values in the `logcounts` assay, and all other cells are assigned `NA` values.
However, as described in the {ref}`ADT data processing section <processing_information:Processed ADT data>`, normalization may fail under certain circumstances.
If a given ADT library failed normalization, it will also have `NA` values in the merged object.

The following additional per-cell data columns for the ADT data can be found in the main `colData` data frame (accessed with `colData(merged_sce)` [as above](#singlecellexperiment-cell-metrics)).

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `altexps_adt_sum`      | UMI count for CITE-seq ADTs                       |
| `altexps_adt_detected` | Number of ADTs detected per cell (ADT count > 0 ) |
| `altexps_adt_percent`  | Percent of `total` UMI count from ADT reads       |


In addition, the following QC statistics from [`DropletUtils::cleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html), which were calculated on individual libraries before objects were merged, can be found in the `colData` of the `"adt"` alternative experiment, accessed with `colData(altExp(merged_sce, "adt"))`.
As in the [main experiment's `colData` slot](#singlecellexperiment-cell-metrics), the `"adt"` alternative experiment `colData` also contains the columns `library_id` and `cell_id`.

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `zero.ambient`   | Indicates whether the cell has zero ambient contamination   |
| `high.controls`  | Indicates whether the cell has unusually high total control counts. Only present if negative/isotype control ADTs were used |
| `sum.controls` |  The sum of counts for all control features. Only present if negative/isotype control ADTs were used |
| `high.ambient`  | Indicates whether the cell has unusually high contamination. Only present if negative/isotype control ADTs were not used |
| `ambient.scale` |  The relative amount of ambient contamination. Only present if negative/isotype control ADTs were not used |
| `discard`  | Indicates whether the cell should be discarded based on ADT QC statistics. The `TRUE` and `FALSE` values in this column correspond, respectively, to values `"Discard"` and `"Keep"` in the `colData(merged_sce)$adt_scpca_filter` column |


Metrics for each of the ADTs assayed can be found as a `DataFrame` stored as `rowData` within the alternative experiment:

```r
rowData(altExp(merged_sce, "adt")) # adt metrics
```

This data frame contains the following columns with statistics for each ADT.
The columns `SCPCL000000-mean` and `SCPCL000000-detected` are present for each library in the merged object that has associated CITE-seq data.

| Column name            | Contents                                                                                                                                                                                                                                 |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `adt_id`               | Name or ID of the ADT                                                                                                                                                                                                                    |
| `target_type`          | Whether each ADT is a target (`target`), negative/isotype control (`neg_control`), or positive control (`pos_control`). If this information was not provided, all ADTs will have been considered targets and will be labeled as `target` |
| `SCPCL000000-mean`     | Mean ADT count across all cells/droplets. Only present for libraries with CITE-seq data                                                                                                                                                  |
| `SCPCL000000-detected` | Percent of cells in which the ADT was detected (ADT count > 0 ). Only present for libraries with CITE-seq data                                                                                                                           |

Finally, additional metadata for ADT processing can be found in the metadata slot of the alternative experiment.

Metadata associated with {ref}`data processing <processing_information:Processing information>` is included in the `metadata` slot as a list.
Similar to the [parent experiment metadata](#singlecellexperiment-experiment-metadata), the `metadata` will contain three items: `library_id`, `sample_id`, and `library_metadata`.
The `library_metadata` contains the same list of metadata as its corresponding parent experiment library's metadata and includes one additional field `ambient_profile`, which holds a named list of the ambient concentrations of each ADT in the given library.

```r
metadata(altExp(merged_sce, "adt")) # adt metadata
```

### Additional SingleCellExperiment components for multiplexed libraries

Merged objects are not available for any projects that contain multiplexed libraries.
This is because there is no guarantee that a unique HTO was used for each sample in a given project, so it would not necessarily be possible to determine which HTO corresponds to which sample in a merged object.

## Components of an AnnData merged object

Before getting started, we highly encourage you to familiarize yourself with the general `AnnData` object structure and functions available as part of the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html).
For the most part, the `AnnData` objects that we provide are formatted to match the expected data format for [`CELLxGENE`](https://cellxgene.cziscience.com/) following [schema version `3.0.0`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md).

To begin, you will need to load the `AnnData` package and read the H5AD file:

```python
import anndata
merged_adata_object = anndata.read_h5ad("SCPCP000000_merged_rna.h5ad")
```

### AnnData expression counts

Merged `AnnData` objects contain two data matrices, each containing RNA-seq expression data for all libraries in the given ScPCA project combined into a single matrix.
The data matrix `raw.X` of the merged `AnnData` object contains the RNA-seq expression data as primary integer counts, and the data matrix `X` contains the RNA-seq expression data as normalized counts.
The data is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a gene.
The `raw.X` and `X` matrices can be accessed with the following python code:

```python
merged_adata_object.raw.X # raw count matrix
merged_adata_object.X # normalized count matrix
```

Column names are cell barcode sequences prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`, and row names are Ensembl gene IDs.
These names can be accessed as with the following python code:


```python
merged_adata_object.obs_names # matrix column names
merged_adata_object.var_names # matrix row names
```

### AnnData cell metrics

Cell metrics calculated from the RNA-seq expression data, which were calculated separately for each library, are stored as a `pandas.DataFrame` in the `.obs` slot.
The slot's row names are cell barcode sequences prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`.

```python
merged_adata_object.obs # cell metrics and metadata
```

All of the per-cell data columns included in the `colData` of the `SingleCellExperiment` merged objects are present in the `.obs` slot of the `AnnData` object.
To see a full description of the included columns, see the sections on [cell metrics](#singlecellexperiment-cell-metrics) and [sample metadata](#singlecellexperiment-sample-metadata) in [`Components of a SingleCellExperiment merged object`](#components-of-a-singlecellexperiment-merged-object).


### AnnData gene information and metrics

Gene information and metrics from the RNA-seq expression data, which were calculated separately for each library, are stored as a `pandas.DataFrame` in the `.var` slot, with the Ensembl ID as the names of the rows.

```python
merged_adata_object.var # gene metrics
```

All of the per-gene data columns included in the `rowData` of the `SingleCellExperiment` objects are present in the `.var` slot of the `AnnData` object.
Note that the `SingleCellExperiment` columns named `SCPCL000000-mean` and `SCPCL000000-detected` are instead named `SCPCL000000.mean` and `SCPCL000000.detected`, respectively, in the merged `AnnData` object.
To see a full description of the included columns, see the [section on gene metrics in `Components of a SingleCellExperiment merged object`](#singlecellexperiment-gene-information-and-metrics).


### AnnData experiment metadata

A partial set of the metadata associated with {ref}`data processing <processing_information:Processing information>` is included in the `.uns` slot of the `AnnData` object as a list.

```python
merged_adata_object.uns # experiment metadata
```

The following items are available in the `.uns` slot:

| Item name                      | Contents                                                                                                                                                                                                          |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `library_id`                   | A list of library IDs which are included in the merged object, each in the form `SCPCL000000`                                                                                                                   |
| `sample_id`                    | A list of sample IDs which are included in the merged object, each in the form `SCPCS000000`                                                                                                                |
| `merged_highly_variable_genes` | A list of highly variable genes used for performing dimensionality reduction on the merged object, determined using `scran::modelGeneVar`, specifying each library as a separate block, and `scran::getTopHVGs` |

Additional experiment metadata is available in the {ref}`metadata TSV file included in the ScPCA Portal download folder <download_files:Metadata>`.

### AnnData dimensionality reduction results

The merged `AnnData` object contains a slot `.obsm` with both principal component analysis (`X_PCA`) and UMAP (`X_UMAP`) results.

For information on how PCA and UMAP results were calculated see the {ref}`section on processed gene expression data <processing_information:Processed gene expression data>`.

The following command can be used to access the PCA and UMAP results:

```python
merged_adata_object.obsm["X_PCA"] # pca results
merged_adata_object.obsm["X_UMAP"] # umap results
```


### Additional AnnData components for CITE-seq libraries (with ADT tags)

ADT data from CITE-seq experiments, when present, is available as a separate `AnnData` object in an H5AD file named with the `_adt.h5ad` suffix.

Merged `AnnData` objects contain two data matrices, each containing CITE-seq expression data for all libraries in a given ScPCA project combined into a single matrix.
The data matrix `raw.X` of the merged `AnnData` object contains the CITE-seq expression data as primary integer counts, and the data matrix `X` contains the RNA-seq expression data as normalized counts.
Note that only cells which are denoted as `"Keep"` in  the `merged_adata_object.uns["adt_scpca_filter"]` column (as described [above](#singlecellexperiment-cell-metrics)) have normalized expression values in the `X` matrix, and all other cells are assigned `NA` values.
The data is stored as a sparse matrix, where each column represents a cell or droplet, and each row represents a single ADT.
The `raw.X` and `X` matrices can be accessed with the following python code:

```python
merged_citeseq_adata_object.raw.X # raw count matrix
merged_citeseq_adata_object.X # normalized count matrix
```

Column names are cell barcode sequences prefixed with the originating library ID, e.g. `SCPCL000000-{barcode}`, and row names are the ADT IDs.

```python
merged_citeseq_adata_object.obs_names # matrix column names
merged_citeseq_adata_object.var_names # matrix row names
```

All of the per-cell data columns included in the `colData` of the `"adt"` alternative experiment in `SingleCellExperiment` merged objects are present in the `.obs` slot of the CITE-seq `AnnData` object.
To see a full description of the included columns, see the section [on additional `SingleCellExperiment` components for CITE-seq libraries](#additional-singlecellexperiment-components-for-cite-seq-libraries-with-adt-tags).


In addition, all of the per-ADT data columns included in the `rowData` of the `"adt"` alternative experiment in `SingleCellExperiment` merged objects are present in the `.var` slot of the CITE-seq `AnnData` object.
To see a full description of the included columns, see the section [on additional `SingleCellExperiment` components for CITE-seq libraries](#additional-singlecellexperiment-components-for-cite-seq-libraries-with-adt-tags).


Finally, the `.uns` slot of the CITE-seq `AnnData` object contains a limited set of experiment metadata information, including a list of library IDs (`library_id`) and sample IDs (`sample_id`) included in the merged object.
