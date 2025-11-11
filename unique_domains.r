#!/usr/bin/env Rscript

# Get command-line argument
N <- commandArgs(trailingOnly = TRUE)[1]

# Construct file name dynamically
input_file <- paste0(N, "_assembled_id_sig_domains")
output_file <- paste0(N, "_assembled_id_sig_unique_domains")

# Read input file
sigd <- read.csv(input_file, header = FALSE, sep = '\t')

# Assign column names
names(sigd) <- c("domain", "transcript_ID", "score", "sspace_boundaries", "cath_resolve_start", "cath_resolve_stop", "cond_evalue", "indp_evalue")

# Ensure domain names are unique
sigd$domains_unique <- make.unique(sigd$domain)

# Select relevant columns
sigd <- sigd[, c("domains_unique", "transcript_ID", "score", "sspace_boundaries", "cath_resolve_start", "cath_resolve_stop", "cond_evalue", "indp_evalue")]

# Write to output file with tab separation
write.table(sigd, output_file, row.names = FALSE, sep = '\t', quote = FALSE)

