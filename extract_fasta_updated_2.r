# Set paths
data_dir <- "./"
output_dir <- "bed_files/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Domain-specific minimum lengths
domain_min_lengths <- c(
  "CIDRa" = 750, "CIDRb" = 765, "CIDRd" = 768, "CIDRg" = 741, "CIDRpam" = 687,
  "DBLa0" = 1225, "DBLa1" = 1167, "DBLa2" = 1252, "DBLb" = 1360, "DBLd" = 1131, "DBLe" = 861, "DBLepam" = 920,
  "DBLg" = 1048, "DBLpam" = 1083, "DBLz" = 1268,
  "NTSA" = 204, "NTSB" = 177, "NTSpam" = 153
)

# Function to extract domain subtype
extract_domain_subtype <- function(domain_name) {
  # Match CIDR, DBL, or NTS subtype anywhere after underscore
  m <- regmatches(domain_name, regexpr("(CIDR[a-zA-Z0-9]+|DBL[a-zA-Z0-9]+|NTS[A-Za-z0-9]+)", domain_name))
  if(length(m) == 0) return(NA)
  sub("", "", m) # remove the leading underscore
}

# Clean subtype for min length lookup (strip numeric suffix)
clean_subtype <- function(subtype) {
  if(is.na(subtype)) return(NA)
  
  # Remove trailing _numbers first (like DBLa0_1 -> DBLa0)
  subtype <- sub("_[0-9]+$", "", subtype)
  
  # If exact match exists in lookup table, keep it
  if(subtype %in% names(domain_min_lengths)) {
    return(subtype)
  }
  
  # Otherwise, remove trailing numbers to get general type
  cleaned <- sub("[0-9]+$", "", subtype)
  return(cleaned)
}

# List all files
file_list <- list.files(path = data_dir, pattern = "_assembled_id_sig_unique_domains$", full.names = TRUE)
cat("Found", length(file_list), "files to process\n")

for(file in file_list) {
  cat("\nProcessing file:", file, "\n")
  if (!file.exists(file) || file.size(file) == 0) {
    cat("File missing or empty, skipping\n")
    next
  }

  # Try reading file
  df <- tryCatch(
    read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      cat("Error reading file with tab, trying space...\n")
      read.delim(file, sep = " ", header = TRUE, stringsAsFactors = FALSE)
    }
  )

  sample_name <- sub("_assembled_id_sig_unique_domains$", "", basename(file))
  cat("Sample name:", sample_name, "\n")

  # Check columns
  if(ncol(df) < 6) {
    cat("Not enough columns, skipping\n")
    next
  }

  # Calculate domain size
  df$domain_size <- as.numeric(df$cath_resolve_stop) - as.numeric(df$cath_resolve_start) + 1 
  # Extract subtype
  df$domain_subtype <- vapply(df$domains_unique, extract_domain_subtype, character(1))
  df$domain_subtype_clean <- vapply(df$domain_subtype, clean_subtype, character(1))

  # Get min lengths
  df$min_length <- domain_min_lengths[df$domain_subtype_clean]

  # Debug info
  cat("Domain subtypes found:", paste(unique(df$domain_subtype), collapse=", "), "\n")
  cat("Sizes range:", min(df$domain_size), "-", max(df$domain_size), "bp\n")

  # Filter by min length
  filtered_df <- df[!is.na(df$min_length) & df$domain_size >= df$min_length, ]
  cat("Passing domains:", paste(unique(filtered_df$domain_subtype), collapse=", "), "\n")

  # Prepare BED columns: transcript_ID, start, stop, domain_name, domain_size
  if(nrow(filtered_df) > 0) {
    bed_df <- filtered_df[, c("transcript_ID", "cath_resolve_start", "cath_resolve_stop", "domains_unique", "domain_size")]
    bed_file <- file.path(output_dir, paste0(sample_name, ".bed"))
    write.table(bed_df, file = bed_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    cat("Wrote BED file:", bed_file, "\n")
  }
}
