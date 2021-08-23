#!/usr/bin/env Rscript
#
# Run spell check and save results
# Adapted from: https://github.com/AlexsLemonade/refinebio-examples/blob/33cdeff66d57f9fe8ee4fcb5156aea4ac2dce07f/scripts/spell-check.R 
# and https://github.com/AlexsLemonade/training-modules/blob/04bea3d2707975e04b57b714ba8b709c77594706/scripts/spell-check.R

library(magrittr)

# Find .git root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Read in dictionary
dictionary <- readLines(file.path(root_dir, 'components', 'dictionary.txt'))

# The only files we want to check are Markdown files
files <- list.files(root_dir, pattern = '\\.md$', recursive = TRUE, full.names = TRUE)


# Run spell check
spelling_errors <- spelling::spell_check_files(files, ignore = dictionary) %>%
  data.frame() %>%
  tidyr::unnest(cols = found) %>%
  tidyr::separate(found, into = c("file", "lines"), sep = ":")

# Print out how many spell check errors
write(nrow(spelling_errors), stdout())

# Save spell errors to file temporarily
readr::write_tsv(spelling_errors, 'spell_check_errors.tsv')
