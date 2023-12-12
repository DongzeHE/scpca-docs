# Single-cell gene expression file contents

Single-cell or single-nuclei gene expression data (unfiltered, filtered, or processed) is provided in two formats:
  - As an RDS file containing a [`SingleCellExperiment` object](http://bioconductor.org/books/3.17/OSCA.intro/the-singlecellexperiment-class.html) for use in R.
  - An HDF5 file containing an [`AnnData` object](https://anndata.readthedocs.io/en/latest/index.html) for use in Python.

These objects contain the expression data, cell and gene metrics, associated metadata, and, in the case of multimodal data like ADTs from CITE-seq experiments, data from additional cell-based assays.
For `SingleCellExperiment` objects, the ADT data will be included as an alternative experiment in the same object containing the primary RNA data.
For `AnnData` objects, the ADT data will be available as a separate object stored in a separate file.
Note that multiplexed sample libraries are only available as `SingleCellExperiment` objects, and are not currently available as `AnnData` objects.

Below we present some details about the specific contents of the objects we provide.

## Components of a SingleCellExperiment object

Before getting started, we highly encourage you to familiarize yourself with the general `SingleCellExperiment` object structure and functions available as part of the [`SingleCellExperiment` package](https://bioconductor.org/packages/3.17/bioc/html/SingleCellExperiment.html) from Bioconductor.

To begin, you will need to load the `SingleCellExperiment` package and read the RDS file:

```r
library(SingleCellExperiment)
sce <- readRDS("SCPCL000000_processed.rds")
```

### SingleCellExperiment expression counts

The `counts` assay of the `SingleCellExperiment` object for single-cell and single-nuclei experiments (for all provided file types) contains the primary RNA-seq expression data as integer counts.
The counts here include reads aligned to both spliced and unspliced cDNA (see the section on {ref}`Post Alevin-fry processing <processing_information:post alevin-fry processing>`).
The data is stored as a sparse matrix, and each column represents a cell or droplet, each row a gene.
Column names are cell barcode sequences and row names are Ensembl gene IDs.
The `counts` assay can be accessed with the following R code:

```r
counts(sce) # counts matrix
```

Additionally, the `spliced` assay contains a counts matrix that includes reads from spliced cDNA only.

### SingleCellExperiment Cell metrics

Cell metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `colData` slot, with the cell barcodes as the names of the rows.

```r
colData(sce) # cell metrics
```

The following per-cell data columns are included for each cell, calculated using the [`scuttle::addPerCellQCMetrics()`](https://rdrr.io/github/LTLA/scuttle/man/addPerCellQCMetrics.html) function.

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sum`                   | UMI count for RNA-seq data                                                                                                                                                                    |
| `detected`              | Number of genes detected (gene count > 0 )                                                                                                                                                    |
| `subsets_mito_sum`      | UMI count of mitochondrial genes                                                                                                                                                              |
| `subsets_mito_detected` | Number of mitochondrial genes detected                                                                                                                                                        |
| `subsets_mito_percent`  | Percent of all UMI counts assigned to mitochondrial genes                                                                                                                                     |
| `total`                 | Total UMI count for RNA-seq data and any alternative experiments (i.e., ADT data from CITE-seq)                                                                                                             |

The following additional per-cell data columns are included in both the `filtered` and `processed` objects.
These columns include metrics calculated by [`miQC`](https://bioconductor.org/packages/release/bioc/html/miQC.html), a package that jointly models proportion of reads belonging to mitochondrial genes and number of unique genes detected to predict low-quality cells.
We also include the filtering results used for the creation of the `processed` objects.
See the description of the {ref}`processed gene expression data <processing_information:Processed gene expression data>` for more information on filtering performed to create the `processed` objects.

| Column name             | Contents                                                                                                                                                                                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `prob_compromised`      | Probability that a cell is compromised (i.e., dead or damaged), as calculated by `miQC`                                                                                                       |
| `miQC_pass`             | Indicates whether the cell passed the default miQC filtering. `TRUE` is assigned to cells with a low probability of being compromised (`prob_compromised` < 0.75) or [sufficiently low mitochondrial content](https://bioconductor.org/packages/release/bioc/vignettes/miQC/inst/doc/miQC.html#preventing-exclusion-of-low-mito-cells)  |
| `scpca_filter` | Labels cells as either `Keep` or `Remove` based on filtering criteria (`prob_compromised` < 0.75 and number of unique genes detected > 200) |
| `adt_scpca_filter` | If CITE-seq was performed, labels cells as either `Keep` or `Remove` based on ADT filtering criteria (`discard = TRUE` as determined by [`DropletUtils::CleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html)) |
| `submitter_celltype_annotation` | If available, cell type annotations obtained from the group that submitted the original data. Cells that the submitter did not annotate are labeled as `NA` |


The `processed` object has one additional `colData` column reflecting cluster assignments.
Further, if cell type annotation was performed, there will be additional columns representing annotation results in the `processed` object's `colData`, as described in the {ref}`cell type annotation processing section <processing_information:Cell type annotation>`.

| Column name             | Contents                                              |
| ----------------------- | ----------------------------------------------------- |
| `cluster`  | Cell cluster identity identified by graph-based clustering |
| `singler_celltype_annotation`  | If cell typing with `SingleR` was performed, the annotated cell type. Cells labeled as `NA` are those which `SingleR` could not confidently annotate |
| `singler_celltype_ontology`  | If cell typing with `SingleR` was performed with ontology labels, the annotated cell type's ontology ID. Cells labeled as `NA` are those which `SingleR` could not confidently annotate |
| `cellassign_celltype_annotation`  | If cell typing with `CellAssign` was performed, the annotated cell type. Cells labeled as `"other"` are those which `CellAssign` could not confidently annotate  |
| `cellassign_max_prediction`  | If cell typing with `CellAssign` was performed, the annotation's prediction score (probability)  |

### SingleCellExperiment gene information and metrics

Gene information and metrics calculated from the RNA-seq expression data are stored as a `DataFrame` in the `rowData` slot, with the Ensembl ID as the names of the rows.

```r
rowData(sce) # gene metrics
```

The following columns are included for all genes.
Metrics were calculated using the [`scuttle::addPerFeatureQCMetrics`](https://rdrr.io/github/LTLA/scuttle/man/addPerFeatureQCMetrics.html) function.

| Column name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `gene_symbol` | [HUGO](https://www.genenames.org) gene symbol, if defined        |
| `mean`        | Mean count across all cells/droplets                             |
| `detected`    | Percent of cells in which the gene was detected (gene count > 0 ) |

### SingleCellExperiment experiment metadata

Metadata associated with {ref}`data processing <processing_information:Processing information>` is included in the `metadata` slot as a list.

```r
metadata(sce) # experiment metadata
```

| Item name           | Contents                                                                                                                       |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `salmon_version`    | Version of `salmon` used for initial mapping                                                                                   |
| `reference_index`   | Transcriptome reference file used for mapping                                                                                  |
| `total_reads`       | Total number of reads processed by `salmon`                                                                                    |
| `mapped_reads`      | Number of reads successfully mapped                                                                                            |
| `mapping_tool`      | Pipeline used for mapping and quantification (`alevin-fry` for all current data in ScPCA)                                      |
| `alevinfry_version` | Version of `alevin-fry` used for mapping and quantification                                                                    |
| `af_permit_type`    | `alevin-fry generate-permit-list` method used for filtering cell barcodes                                                      |
| `af_resolution`     | `alevin-fry quant` resolution mode used                                                                                        |
| `usa_mode`          | Boolean indicating whether quantification was done using `alevin-fry` USA mode                                                 |
| `af_num_cells`      | Number of cells reported by `alevin-fry`                                                                                       |
| `tech_version`      | A string indicating the technology and version used for the single-cell library, such as 10Xv2, 10Xv3, or 10Xv3.1              |
| `assay_ontology_term_id` | A string indicating the [Experimental Factor Ontology](https://www.ebi.ac.uk/ols/ontologies/efo) term id associated with the `tech_version`  |
| `seq_unit`         | `cell` for single-cell samples or `nucleus` for single-nucleus samples                                                          |
| `transcript_type`   | Transcripts included in gene counts: `spliced` for single-cell samples and `unspliced` for single-nuclei                       |
| `sample_metadata`   | Data frame containing metadata for each sample included in the library (see the [`Sample metadata` section below](#singlecellexperiment-sample-metadata)) |
| `miQC_model`        | The model object that `miQC` fit to the data and was used to calculate `prob_compromised`. Only present for `filtered` objects |
| `filtering_method`  | The method used for cell filtering. One of `emptyDrops`, `emptyDropsCellRanger`, or `UMI cutoff`. Only present for `filtered` objects |
| `umi_cutoff`        | The minimum UMI count per cell used as a threshold for removing empty droplets. Only present for `filtered` objects where the `filtering_method` is `UMI cutoff` |
| `prob_compromised_cutoff`        | The minimum cutoff for the probability of a cell being compromised, as calculated by `miQC`. Only present for `filtered` objects |
| `scpca_filter_method`        | Method used by the Data Lab to filter low quality cells prior to normalization. Either `miQC` or `Minimum_gene_cutoff`  |
| `adt_scpca_filter_method`        | If CITE-seq was performed, the method used by the Data Lab to identify cells to be filtered prior to normalization, based on ADT counts. Either `cleanTagCounts with isotype controls` or `cleanTagCounts without isotype controls`. If filtering failed (i.e. `DropletUtils::cleanTagCounts()` could not reliably determine which cells to filter), the value will be `No filter` |
| `min_gene_cutoff`        | The minimum cutoff for the number of unique genes detected per cell. Only present for `filtered` objects |
| `normalization`        | The method used for normalization of raw RNA counts. Either `deconvolution`, described in [Lun, Bach, and Marioni (2016)](https://doi.org/10.1186/s13059-016-0947-7), or `log-normalization`. Only present for `processed` objects |
| `adt_normalization`        | If CITE-seq was performed, the method used for normalization of raw ADT counts. Either `median-based` or  `log-normalization`, as explained in the {ref}`processed ADT data section <processing_information:Processed ADT data>`. Only present for `processed` objects |
| `highly_variable_genes`        | A list of highly variable genes used for dimensionality reduction, determined using `scran::modelGeneVar` and `scran::getTopHVGs`. Only present for `processed` objects |
| `cluster_algorithm` | The algorithm used to perform graph-based clustering of cells. Only present for `processed` objects |
| `cluster_weighting` | The weighting approach used during graph-based clustering. Only present for `processed` objects |
| `cluster_nn`        | The nearest neighbor parameter value used for the graph-based clustering. Only present for `processed` objects |
| `celltype_methods` | If cell type annotation was performed, a vector of the methods used for annotation. May include `"submitter"`, `"singler"` and/or `"cellassign"`. If submitter cell-type annotations are available, this metadata item will be present in all objects. Otherwise, this item will only be in `processed` objects |
| `singler_results` | If cell typing with `SingleR` was performed, the full result object returned by `SingleR` annotation. Only present for `processed` objects |
| `singler_reference` | If cell typing with `SingleR` was performed, the name of the reference dataset used for annotation. Only present for `processed` objects |
| `singler_reference_label` | If cell typing with `SingleR` was performed, the name of the label in the reference dataset used for annotation. Only present for `processed` objects |
| `singler_reference_source`  | If cell typing with `SingleR` was performed, the source of the reference dataset (default is [`celldex`](http://bioconductor.org/packages/release/data/experiment/html/celldex.html)). Only present for `processed` objects |
| `singler_reference_version`  | If cell typing with `SingleR` was performed, the version of `celldex` used to create the reference dataset source, with periods replaced as dashes (`-`). Only present for `processed` objects |
| `cellassign_predictions` | If cell typing with `CellAssign` was performed, the full matrix of predictions across cells and cell types. Only present for `processed` objects |
| `cellassign_reference` | If cell typing with `CellAssign` was performed, the reference name as established by the Data Lab used for cell type annotation. Only present for `processed` objects |
| `cellassign_reference_organs` | If cell typing with `CellAssign` was performed, a comma-separated list of organs and/or tissue groupings for which marker genes were obtained to create the reference. Only present for `processed` objects |
| `cellassign_reference_source`  | If cell typing with `CellAssign` was performed, the source of the reference dataset (default is [`PanglaoDB`](https://panglaodb.se/)). Only present for `processed` objects |
| `cellassign_reference_version`  | If cell typing with `CellAssign` was performed, the version of the reference dataset source. For references obtained from `PanglaoDB`, the version scheme is a date in ISO8601 format. Only present for `processed` objects |


### SingleCellExperiment sample metadata

Relevant sample metadata is available as a data frame stored in the `metadata(sce)$sample_metadata` slot of the `SingleCellExperiment` object.
Each row in the data frame will correspond to a sample present in the library.
The following columns are included in the sample metadata data frame for all libraries.

| Column name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `sample_id`   | Sample ID in the form `SCPCS000000`                            |
| `library_id`   | Library ID in the form `SCPCL000000`                             |
| `particpant_id`  | Unique id corresponding to the donor from which the sample was obtained |
| `submitter_id`    | Original sample identifier from submitter                      |
| `submitter`       | Submitter name/id                                              |
| `age`             | Age at time sample was obtained                                |
| `sex`             | Sex of patient that the sample was obtained from               |
| `diagnosis`       | Tumor type                                                     |
| `subdiagnosis`    | Subcategory of diagnosis or mutation status (if applicable)    |
| `tissue_location` | Where in the body the tumor sample was located                 |
| `disease_timing`  | At what stage of disease the sample was obtained, either diagnosis or recurrence |
| `organism`         | The organism the sample was obtained from (e.g., `Homo_sapiens`) |
| `development_stage_ontology_term_id` | [`HsapDv` ontology](http://obofoundry.org/ontology/hsapdv.html) term indicating developmental stage. If unavailable, `unknown` is used  |
| `sex_ontology_term_id` | [`PATO`](http://obofoundry.org/ontology/pato.html) term referring to the sex of the sample. If unavailable, `unknown` is used |
| `organism_ontology_id` | [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) term for organism, e.g. [`NCBITaxon:9606`](http://purl.obolibrary.org/obo/NCBITaxon_9606) |
| `self_reported_ethnicity_ontology_term_id` | For _Homo sapiens_, a [`Hancestro` term](http://obofoundry.org/ontology/hancestro.html). `multiethnic` indicates more than one ethnicity is reported. `unknown` indicates unavailable ethnicity and `NA` is used for all other organisms  |
| `disease_ontology_term_id` | [`MONDO`](http://obofoundry.org/ontology/mondo.html) term indicating disease type. [`PATO:0000461`](http://purl.obolibrary.org/obo/PATO_0000461) indicates normal or healthy tissue. If unavailable, `NA` is used  |
| `tissue_ontology_term_id`| [`UBERON`](http://obofoundry.org/ontology/uberon.html) term indicating tissue of origin. If unavailable, `NA` is used |

For some libraries, the sample metadata may also include additional metadata specific to the disease type and experimental design of the project.
Examples of this include treatment or outcome.

### SingleCellExperiment dimensionality reduction results

In the RDS file containing the processed `SingleCellExperiment` object only (`_processed.rds`), the `reducedDim` slot of the object will be occupied with both principal component analysis (`PCA`) and `UMAP` results.
For all other files, the `reducedDim` slot will be empty as no dimensionality reduction was performed.

PCA results were calculated using `scater::runPCA()`, using only highly variable genes.
The list of highly variable genes used was selected using `scran::modelGeneVar` and `scran::getTopHVGs`, and are stored in the `SingleCellExperiment` object in `metadata(sce)$highly_variable_genes`.
The following command can be used to access the PCA results:

```r
reducedDim(sce, "PCA")
```

UMAP results were calculated using `scater::runUMAP()`, with the PCA results as input rather than the full gene expression matrix.
The following command can be used to access the UMAP results:

```r
reducedDim(sce, "UMAP")
```

### Additional SingleCellExperiment components for CITE-seq libraries (with ADT tags)

ADT data from CITE-seq experiments, when present, is included within the `SingleCellExperiment` as an "Alternative Experiment" named `"adt"` , which can be accessed with the following command:

```r
altExp(sce, "adt") # adt experiment
```

Within this, the main expression matrix is again found in the `counts` assay and the normalized expression matrix is found in the `logcounts` assay.
For each assay, each column corresponds to a cell or droplet (in the same order as the parent `SingleCellExperiment`) and each row corresponds to an antibody derived tag (ADT).
Column names are again cell barcode sequences and row names are the antibody targets for each ADT.

Note that only cells which are denoted as "Keep" in  the `colData(sce)$adt_scpca_filter` column (as described [above](#singlecellexperiment-cell-metrics)) have normalized expression values in the `logcounts` assay, and all other cells are assigned `NA` values.
However, as described in the {ref}`processed ADT data section <processing_information:Processed ADT data>`, normalization may fail under certain circumstances, in which case there will be no `logcounts` normalized expression matrix present in the alternative experiment.

The following additional per-cell data columns for the ADT data can be found in the main `colData` data frame (accessed with `colData(sce)` [as above](#singlecellexperiment-cell-metrics)).

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `altexps_adt_sum`      | UMI count for CITE-seq ADTs                       |
| `altexps_adt_detected` | Number of ADTs detected per cell (ADT count > 0 ) |
| `altexps_adt_percent`  | Percent of `total` UMI count from ADT reads       |


In addition, the following QC statistics from [`DropletUtils::cleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html) can be found in the `colData` of the `"adt"` alternative experiment, accessed with `colData(altExp(sce, "adt"))`.

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `zero.ambient`   | Indicates whether the cell has zero ambient contamination   |
| `sum.controls` |  The sum of counts for all control features. Only present if negative/isotype control ADTs are present |
| `high.controls`  | Indicates whether the cell has unusually high total control counts. Only present if negative/isotype control ADTs are present |
| `ambient.scale` |  The relative amount of ambient contamination. Only present if negative control ADTs are _not_ present |
| `high.ambient`  | Indicates whether the cell has unusually high contamination. Only present if negative/isotype control ADTs are _not_ present |
| `discard`  | Indicates whether the cell should be discarded based on QC statistics |


Metrics for each of the ADTs assayed can be found as a `DataFrame` stored as `rowData` within the alternative experiment:

```r
rowData(altExp(sce, "adt")) # adt metrics
```

This data frame contains the following columns with statistics for each ADT:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `mean`      | Mean ADT count across all cells/droplets                       |
| `detected`  | Percent of cells in which the ADT was detected (ADT count > 0 ) |
| `target_type` | Whether each ADT is a target (`target`), negative/isotype control (`neg_control`), or positive control (`pos_control`). If this information was not provided, all ADTs will have been considered targets and will be labeled as `target` |

Finally, additional metadata for ADT processing can be found in the metadata slot of the alternative experiment.
This metadata slot has the same contents as the [parent experiment metadata](#singlecellexperiment-experiment-metadata), along with one additional field, `ambient_profile`, which holds a list of the ambient concentrations of each ADT.

```r
metadata(altExp(sce, "adt")) # adt metadata
```

### Additional SingleCellExperiment components for multiplexed libraries

Multiplexed libraries will contain a number of additional components and fields.

Hashtag oligo (HTO) quantification for each cell is included within the `SingleCellExperiment` as an "Alternative Experiment" named `"cellhash"` , which can be accessed with the following command:

```r
altExp(sce, "cellhash") # hto experiment
```

Within this, the main data matrix is again found in the `counts` assay, with each column corresponding to a cell or droplet (in the same order as the parent `SingleCellExperiment`) and each row corresponding to a hashtag oligo (HTO).
Column names are again cell barcode sequences and row names the HTO ids for all assayed HTOs.

The following additional per-cell data columns for the cellhash data can be found in the main `colData` data frame (accessed with `colData(sce)` [as above](#singlecellexperiment-cell-metrics)).

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `altexps_cellhash_sum`    | UMI count for cellhash HTOs                       |
| `altexps_cellhash_detected` | Number of HTOs detected per cell (HTO count > 0 ) |
| `altexps_cellhash_percent`  | Percent of `total` UMI count from HTO reads       |

Metrics for each of the HTOs assayed can be found as a `DataFrame` stored as `rowData` within the alternative experiment:

```r
rowData(altExp(sce, "cellhash")) # hto metrics
```

This data frame contains the following columns with statistics for each HTO:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `mean`      | Mean HTO count across all cells/droplets                       |
| `detected`  | Percent of cells in which the HTO was detected (HTO count > 0 ) |
| `sample_id` | Sample ID for this library that corresponds to the HTO (only present in `_filtered.rds` files) |

Note that in the unfiltered `SingleCellExperiment` objects, this may include hashtag oligos that do not correspond to any included sample, but were part of the reference set used for mapping.

### Demultiplexing results

Demultiplexing results are included only in the `_filtered.rds` files.
The demultiplexing methods applied for these files are as described in the {ref}`multiplex data processing section <processing_information:Multiplexed libraries>`.

Demultiplexing analysis adds the following additional fields to the `colData(sce)` data frame:

| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `hashedDrops_sampleid`  | Most likely sample as called by `DropletUtils::hashedDrops` |
| `HTODemux_sampleid`  | Most likely sample as called by `Seurat::HTODemux` |
| `vireo_sampleid`  | Most likely sample as called by `vireo` (genetic demultiplexing) |

### Additional demultiplexing statistics

Each demultiplexing method generates additional statistics specific to the method that you may wish to access, including probabilities, alternative calls, and potential doublet information.

For methods that rely on the HTO data, these statistics are found in the `colData(altExp(sce, "cellhash"))` data frame;
`DropletUtils::hashedDrops()` statistics have the prefix `hashedDrops_` and `Seurat::HTODemux()` statistics have the prefix `HTODemux`.

Genetic demultiplexing statistics are found in the main `colData(sce)` data frame, with the prefix `vireo_`.

## Components of an AnnData object

Before getting started, we highly encourage you to familiarize yourself with the general `AnnData` object structure and functions available as part of the [`AnnData` package](https://anndata.readthedocs.io/en/latest/index.html).
For the most part, the `AnnData` objects that we provide are formatted to match the expected data format for [`CELLxGENE`](https://cellxgene.cziscience.com/) following [schema version `3.0.0`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md).

To begin, you will need to load the `AnnData` package and read the HDF5 file:

```python
import anndata
adata_object = anndata.read_h5ad("SCPCL000000_processed_rna.hdf5")
```

### AnnData expression counts

The data matrix, `X`, of the `AnnData` object for single-cell and single-nuclei experiments contains the primary RNA-seq expression data as integer counts in both the unfiltered (`_unfiltered_rna.hdf5`) and filtered (`_filtered_rna.hdf5`) objects.
The data is stored as a sparse matrix, and each column represents a cell or droplet, each row a gene.
Column names are cell barcode sequences and row names are Ensembl gene IDs.
The `X` matrix can be accessed with the following python code:

```python
adata_object.X # raw count matrix
```

In processed objects _only_ (`_processed_rna.hdf5`), the data matrix `X` contains the normalized data and the primary data can be found in `raw.X`.
The counts in the processed object can be accessed with the following python code:

```python
adata_object.raw.X # raw count matrix
adata_object.X # normalized count matrix
```

### AnnData cell metrics

Cell metrics calculated from the RNA-seq expression data are stored as a `pandas.DataFrame` in the `.obs` slot, with the cell barcodes as the names of the rows.

```python
adata_object.obs # cell metrics
```

All of the per-cell data columns included in the `colData` of the `SingleCellExperiment` objects are present in the `.obs` slot of the `AnnData` object.
To see a full description of the included columns, see the [section on cell metrics in Components of a SingleCellExperiment object](#singlecellexperiment-cell-metrics).

The `AnnData` object also includes the following additional cell-level metadata columns:

| Column name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `sample_id`   | Sample ID in the form `SCPCS000000`                            |
| `library_id`   | Library ID in the form `SCPCL000000`                             |
| `assay_ontology_term_id` | A string indicating the [Experimental Factor Ontology](https://www.ebi.ac.uk/ols/ontologies/efo) term id associated with the technology and version used for the single-cell library, such as 10Xv2, 10Xv3, or 10Xv3.1 |
| `suspension_type`         | `cell` for single-cell samples or `nucleus` for single-nucleus samples  |
| `particpant_id`  | Unique id corresponding to the donor from which the sample was obtained |
| `submitter_id`    | Original sample identifier from submitter                      |
| `submitter`       | Submitter name/id                                              |
| `age`             | Age at time sample was obtained                                |
| `sex`             | Sex of patient that the sample was obtained from               |
| `diagnosis`       | Tumor type                                                     |
| `subdiagnosis`    | Subcategory of diagnosis or mutation status (if applicable)    |
| `tissue_location` | Where in the body the tumor sample was located                 |
| `disease_timing`  | At what stage of disease the sample was obtained, either diagnosis or recurrence |
| `organism`         | The organism the sample was obtained from (e.g., `Homo_sapiens`) |
| `development_stage_ontology_term_id` | [`HsapDv` ontology](http://obofoundry.org/ontology/hsapdv.html) term indicating developmental stage. If unavailable, `unknown` is used  |
| `sex_ontology_term_id` | [`PATO`](http://obofoundry.org/ontology/pato.html) term referring to the sex of the sample. If unavailable, `unknown` is used |
| `organism_ontology_id` | [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) term for organism, e.g. [`NCBITaxon:9606`](http://purl.obolibrary.org/obo/NCBITaxon_9606) |
| `self_reported_ethnicity_ontology_term_id` | For _Homo sapiens_, a [`HANCESTRO` term](http://obofoundry.org/ontology/hancestro.html). `multiethnic` indicates more than one ethnicity is reported. `unknown` indicates unavailable ethnicity, and `NA` is used for all other organisms  |
| `disease_ontology_term_id` | [`Mondo`](http://obofoundry.org/ontology/mondo.html) term indicating disease type. [`PATO:0000461`](http://purl.obolibrary.org/obo/PATO_0000461) indicates normal or healthy tissue. If unavailable, `NA` is used  |
| `tissue_ontology_term_id`| [`Uberon`](http://obofoundry.org/ontology/uberon.html) term indicating tissue of origin. If unavailable, `NA` is used |
| `is_primary_data` | Set to `FALSE` for all libraries to reflect that all libraries were obtained from external investigators. Required by `CELLxGENE`             |


### AnnData gene information and metrics

Gene information and metrics calculated from the RNA-seq expression data are stored as a `pandas.DataFrame` in the `.var` slot, with the Ensembl ID as the names of the rows.

```python
adata_object.var # gene metrics
```

All of the per-gene data columns included in the `rowData` of the `SingleCellExperiment` objects are present in the `.var` slot of the `AnnData` object.
To see a full description of the included columns, see the [section on gene metrics in `Components of a SingleCellExperiment object`](#singlecellexperiment-gene-information-and-metrics).

The `AnnData` object also includes the following additional gene-level metadata columns:

| Column name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `is_feature_filtered` | Boolean indicating if the gene or feature is filtered out in the normalized matrix but is present in the raw matrix     |


### AnnData experiment metadata

Metadata associated with {ref}`data processing <processing_information:Processing information>` is included in the `.uns` slot as a list.

```python
adata_object.uns # experiment metadata
```

All of the object metadata included in `SingleCellExperiment` objects are present in the `.uns` slot of the `AnnData` object.
To see a full description of the included columns, see the [section on experiment metadata in `Components of a SingleCellExperiment object`](#singlecellexperiment-experiment-metadata).
The only exception is that the `AnnData` object _does not_ contain the `sample_metadata` item in the `.uns` slot.
Instead, the contents of the `sample_metadata` data frame are stored in the cell-level metadata (`.obs`).

The `AnnData` object also includes the following additional items in the `.uns` slot:

| Item name   | Contents                                                         |
| ------------- | ---------------------------------------------------------------- |
| `schema_version` | CZI schema version used for `AnnData` formatting |


### AnnData dimensionality reduction results

The HDF5 file containing the processed `AnnData` object (`_processed_rna.hdf5`) contains a slot `.obsm` with both principal component analysis (`X_PCA`) and UMAP (`X_UMAP`) results.
For all other HDF5 files, the `.obsm` slot will be empty as no dimensionality reduction was performed.

For information on how PCA and UMAP results were calculated see the {ref}`section on processed gene expression data <processing_information:Processed gene expression data>`.

The following command can be used to access the PCA and UMAP results:

```python
adata_object.obsm["X_PCA"] # pca results
adata_object.obsm["X_UMAP"] # umap results
```

### Additional AnnData components for CITE-seq libraries (with ADT tags)

ADT data from CITE-seq experiments, when present, is available as a separate `AnnData` object (HDF5 file).
All files containing ADT data will contain the `_adt.hdf5` suffix.

The data matrix, `X`, of the `AnnData` objects contain the primary ADT expression data as integer counts in both the unfiltered (`_unfiltered_adt.hdf5`) and filtered (`_filtered_adt.hdf5`) objects.
Each column corresponds to a cell or droplet (in the same order as the main `AnnData` object), and each row corresponds to an antibody derived tag (ADT).
Column names are again cell barcode sequences and row names are the antibody targets for each ADT.

As with the RNA `AnnData` objects, in processed objects _only_ (`_processed_adt.hdf5`), the data matrix `X` contains the normalized ADT counts and the primary data can be found in `raw.X`.
Note that only cells which are denoted as `"Keep"` in  the `adata_obj.uns["adt_scpca_filter"]` column (as described [above](#singlecellexperiment-cell-metrics)) have normalized expression values in the `X` matrix, and all other cells are assigned `NA` values.
However, as described in the {ref}`processed ADT data section <processing_information:Processed ADT data>`, normalization may fail under certain circumstances.
In such cases the `AnnData` object will not contain a normalized expression matrix, but the primary data will still be stored in `X`.

In addition, the following QC statistics from [`DropletUtils::cleanTagCounts()`](https://rdrr.io/github/MarioniLab/DropletUtils/man/cleanTagCounts.html) can be found in the `obs` slot of each ADT-specific `AnnData` object.

| Column name                | Contents                                          |
| -------------------------- | ------------------------------------------------- |
| `zero.ambient`   | Indicates whether the cell has zero ambient contamination   |
| `sum.controls` |  The sum of counts for all control features. Only present if negative/isotype control ADTs are present |
| `high.controls`  | Indicates whether the cell has unusually high total control counts. Only present if negative/isotype control ADTs are present |
| `ambient.scale` |  The relative amount of ambient contamination. Only present if negative control ADTs are _not_ present |
| `high.ambient`  | Indicates whether the cell has unusually high contamination. Only present if negative/isotype control ADTs are _not_ present |
| `discard`  | Indicates whether the cell should be discarded based on QC statistics |


Metrics for each of the ADTs assayed can be found as a `pandas.DataFrame` in the `.var` slot within the `AnnData` object:


| Column name | Contents                                                       |
| ----------- | -------------------------------------------------------------- |
| `mean`      | Mean ADT count across all cells/droplets                       |
| `detected`  | Percent of cells in which the ADT was detected (ADT count > 0 ) |
| `target_type` | Whether each ADT is a target (`target`), negative/isotype control (`neg_control`), or positive control (`pos_control`). If this information was not provided, all ADTs will have been considered targets and will be labeled as `target` |

Finally, additional metadata for ADT processing can be found in the `.uns` slot of the `AnnData` object.
This metadata slot has the same contents as the [RNA experiment metadata](#anndata-experiment-metadata), along with one additional field, `ambient_profile`, which holds a list of the ambient concentrations of each ADT.
