library(magrittr)
library(ggplot2)
library(SingleCellExperiment)
library(patchwork)

# path to results files with sces and qc dataframes from running benchmarking_generate_qc_df.R
base_dir <- here::here()
file_dir <- file.path(base_dir, "data", "results")

# sce objects
cr_like_em_file <- file.path(file_dir, "alevin-fry-cr-like-em-emptydrops-200-sces.rds")
cellranger_file <- file.path(file_dir, "cellranger_sces.rds")

# plots to save
umi_boxplot <- file.path("..", "images", "total_umi_per_cell.pdf")
gene_boxplot <- file.path("..","images", "total_genes_per_cell.pdf")

# read in mito gene list
library_data_dir <- file.path(base_dir, 'sample-info')
mito_file <- file.path(library_data_dir, "Homo_sapiens.GRCh38.103.mitogenes.txt")
mito_genes <- readr::read_tsv(mito_file, col_names = "gene_id")
mito_genes <- mito_genes %>%
  dplyr::pull(gene_id) %>%
  unique()

# read in cr_like_em_sce and cellranger sces and add colData QC 
cr_like_em_sce <- readr::read_rds(cr_like_em_file) %>%
  purrr::map(scpcaTools::add_cell_mito_qc, mito = mito_genes)
cellranger_sce <- readr::read_rds(cellranger_file) %>%
  purrr::map(scpcaTools::add_cell_mito_qc, mito = mito_genes)

# make combined list of all sces 
all_sces_list <- list(cr_like_em_sce,cellranger_sce)
names(all_sces_list) <- c("Alevin-fry", "Cell Ranger")

# make list of samples to use for assigning samples to cell or nuclei 
single_cell <- c("SCPCR000126", "SCPCR000127")
single_nuclei <- c("SCPCR000118", "SCPCR00119", "SCPCR000220", "SCPCR000221")

# grab coldata from all samples and create a dataframe used for plotting
coldata_df <- purrr::map_df(
  all_sces_list,
  ~ purrr::map_df(.x, scpcaTools::coldata_to_df, .id = "quant_id"),
  .id = "tool"
) %>% 
  # create new columns with just sample ID and then seq unit 
  dplyr::mutate(sample = stringr::word(quant_id, 1, sep = "-"),
                seq_unit = dplyr::case_when(sample %in% single_cell ~ "Cell",
                                            sample %in% single_nuclei ~ "Nuclei")) %>%
  # remove unwanted samples that were used in benchmarking
  dplyr::filter(!sample %in% c("SCPCR000118", "SCPCR000119"))

# filter for cells that are found in both alevin + cellranger
cell_counts <- coldata_df %>%  
  dplyr::count(cell_id, sample)

common_cells <- cell_counts %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(cell_id)

coldata_common <- coldata_df %>%
  dplyr::filter(cell_id %in% common_cells) 

# create two separate dataframes used for plotting
cell_coldata_common <- coldata_common %>%
  dplyr::filter(seq_unit == "Cell")

nuclei_coldata_common <- coldata_common %>%
  dplyr::filter(seq_unit == "Nuclei")

# create combined UMI per cell plot 
umi_cell <- ggplot(cell_coldata_common, aes(x = tool, y = sum, fill = seq_unit)) + 
  geom_boxplot() + 
  facet_wrap(~ sample) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 8)) +
  labs(x = "",
       y = "Total UMI per cell") +
  scale_fill_manual(values = c("red", "lightblue"))

umi_nuclei <- ggplot(nuclei_coldata_common, aes(x = tool, y = sum, fill = seq_unit)) + 
  geom_boxplot() + 
  facet_wrap(~ sample) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 8)) +
  labs(x = "",
       y = "") + 
  coord_cartesian(ylim = c(0,40000)) +
  scale_fill_manual(values = c("lightblue", "red"))

# combine plots and add combined legend 
combined_umi_plot <- umi_cell + umi_nuclei + plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.title = element_blank())

# write combined plot to pdf
pdf(file = umi_boxplot, width = 7, height = 4)
combined_umi_plot
dev.off()

# do the same thing for genes detected per cell 
genes_cell <- ggplot(cell_coldata_common, aes(x = tool, y = detected, fill = seq_unit)) + 
  geom_boxplot() + 
  facet_wrap(~ sample) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 8)) +
  labs(x = "",
       y = "Total genes detected per cell") +
  scale_fill_manual(values = c("red", "lightblue"))

genes_nuclei <- ggplot(nuclei_coldata_common, aes(x = tool, y = detected, fill = seq_unit)) + 
  geom_boxplot() + 
  facet_wrap(~ sample) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 8)) +
  labs(x = "",
       y = "") + 
  coord_cartesian(ylim = c(0,10000)) +
  scale_fill_manual(values = c("lightblue", "red"))

combined_gene_plot <- genes_cell + genes_nuclei + plot_layout(guides = "collect") & 
  theme(legend.position = "top", legend.title = element_blank())

pdf(file = gene_boxplot, width = 7, height = 4)
combined_gene_plot
dev.off()