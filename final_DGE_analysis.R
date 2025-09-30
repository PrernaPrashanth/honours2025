# Load libraries
library(edgeR); library(limma); library(quadprog); library(ruv)
library(ggplot2); library(reshape2); library(EDASeq); library(RColorBrewer)
library(VennDiagram); library(data.table); library(dplyr); library(stringr)
library(knitr); library(pheatmap); library(ggfortify); library(readxl)
library(ggrepel); library(tidyr); library(writexl); library(gridExtra)

wd = "./"
colors <- c("#0571b0","#ca0020")
setwd("~/Downloads")

# Load and process gene mappings
geneID_mappings_raw <- fread("geneID_mappings.txt", data.table = FALSE, skip = 1)
colnames(geneID_mappings_raw) <- c("GeneID", "GenomicSeqID", "UniProtID", "EntrezID", "Symbol", "PreviousIDs", "V7")
s <- str_split(as.character(geneID_mappings_raw$PreviousIDs), ",")
geneID_mappings <- data.frame(current = rep(geneID_mappings_raw$GeneID, lengths(s)),
                              old = str_trim(unlist(s)))

# Load counts and annotations
counts_raw <- read.table("SM_domain_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$geneid; counts_raw$geneid <- NULL
colnames(counts_raw) <- gsub("_counts$", "", colnames(counts_raw))

annotations <- read.table("annotation_simple.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                          col.names = c("GeneID", "Chr", "Start", "End", "Strand", "Length"))

# Fix duplicate isoforms
pf3d7_rows <- grepl("^PF3D7_", rownames(counts_raw))
base_genes <- gsub("\\.\\d+$", "", rownames(counts_raw))
pf3d7_summed <- aggregate(counts_raw[pf3d7_rows, ], by = list(base_genes[pf3d7_rows]), FUN = sum)
rownames(pf3d7_summed) <- pf3d7_summed$Group.1; pf3d7_summed$Group.1 <- NULL
counts_final <- rbind(pf3d7_summed, counts_raw[!pf3d7_rows, ])

# Create combined annotations
pf3d7_ids <- rownames(counts_final)[grepl("^PF3D7_", rownames(counts_final))]
domain_ids <- rownames(counts_final)[!grepl("^PF3D7_", rownames(counts_final))]

gene_annot <- annotations[annotations$GeneID %in% pf3d7_ids, ]
gene_annot$Feature_Type <- "Core_Gene"; gene_annot$Domain_Family <- NA

domain_annot <- data.frame(GeneID = domain_ids, Chr = "Domain", Start = NA, End = NA, 
                           Strand = NA, Length = NA, Feature_Type = "Protein_Domain",
                           Domain_Family = case_when(grepl("^CIDR", domain_ids) ~ "CIDR",
                                                     grepl("^DBL", domain_ids) ~ "DBL",
                                                     TRUE ~ "Other_Domain"))

combined_annotations <- rbind(gene_annot, domain_annot)
annotations_ordered <- combined_annotations[match(rownames(counts_final), combined_annotations$GeneID), ]

# Create DGEList
x <- DGEList(counts = counts_final, genes = annotations_ordered)

# Fix missing annotations
missing_indices <- which(is.na(x$genes$GeneID))
if(length(missing_indices) > 0) {
  for(i in missing_indices) {
    gene_name <- rownames(x$counts)[i]
    x$genes[i, "GeneID"] <- gene_name
    if(grepl("^PF3D7_", gene_name)) {
      x$genes[i, c("Chr", "Start", "End", "Strand", "Length")] <- NA
      x$genes[i, "Feature_Type"] <- "Core_Gene"; x$genes[i, "Domain_Family"] <- NA
    } else {
      x$genes[i, "Chr"] <- "Domain"; x$genes[i, c("Start", "End", "Strand", "Length")] <- NA
      x$genes[i, "Feature_Type"] <- "Protein_Domain"
      x$genes[i, "Domain_Family"] <- case_when(grepl("^CIDR", gene_name) ~ "CIDR",
                                               grepl("^DBL", gene_name) ~ "DBL",
                                               TRUE ~ "Other_Domain")
    }
  }
}

# Filter and create DGE object
keep <- rowSums(cpm(x) > 2) >= 10
x <- x[keep, , keep.lib.sizes = FALSE]

library(Biostrings)

setwd("~/Desktop/extracted_fasta/")
fasta_files <- list.files(pattern = "*.fasta")

# Create a mapping of sample -> domain -> length
sample_domain_lengths <- list()

for(file in fasta_files) {
  # Extract sample name
  sample_name <- gsub("\\.extracted\\.fasta$", "", basename(file))
  
  # Read sequences
  sequences <- readDNAStringSet(file)
  lengths <- width(sequences)
  seq_names <- names(sequences)
  
  # Extract domain family names from headers
  domain_families <- gsub(".*_scaffold[0-9]+_([^:]+)::.*", "\\1", seq_names)
  
  # Handle duplicate domain names by averaging lengths
  unique_domains <- unique(domain_families)
  domain_lengths <- sapply(unique_domains, function(domain) {
    mean(lengths[domain_families == domain])
  })
  
  sample_domain_lengths[[sample_name]] <- domain_lengths
}

# Go back to original directory
setwd("~/Downloads/")

# Update lengths in DGE object based on sample-specific information
updates_made <- 0
total_possible_updates <- 0

for(j in 1:ncol(x$counts)) {
  sample_name <- colnames(x$counts)[j]
  
  if(sample_name %in% names(sample_domain_lengths)) {
    sample_lengths <- sample_domain_lengths[[sample_name]]
    
    for(i in 1:nrow(x$genes)) {
      gene_id <- x$genes$GeneID[i]
      if(is.na(x$genes$Length[i])) {
        total_possible_updates <- total_possible_updates + 1
        if(gene_id %in% names(sample_lengths)) {
          x$genes$Length[i] <- round(sample_lengths[gene_id])
          updates_made <- updates_made + 1
        }
      }
    }
  }
}

# Check final status
var_genes_na_final <- sum(is.na(x$genes$Length) & !grepl("^PF3D7_", x$genes$GeneID))
cat("Var genes with NA length after update:", var_genes_na_final, "\n")


###################################################################################################################
# READ COUNTS VISUALIZATION
###################################################################################################################

# Calculate total read counts per sample
sample_read_counts <- colSums(x$counts)

# Prepare data for plotting
bplot <- reshape2::melt(sample_read_counts)
colnames(bplot) <- c("Read.Counts")
bplot$Sample <- names(sample_read_counts)

# Plot sample read depths
gg_reads <- ggplot(bplot, aes(x = factor(Sample), y = Read.Counts)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_sqrt(breaks = c(0, 10000, 1000000, 10000000, 100000000)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  labs(x = 'Sample', y = 'Read Count') +
  geom_hline(aes(yintercept = 1e6), col = "red")

print(gg_reads)
ggsave("sample_read_depths.pdf", gg_reads, width = 16, height = 8)

# Print sample filtering summary
cat("Number of samples below 1,000,000 reads:", sum(sample_read_counts < 1e6), "\n")
cat("Samples below threshold:\n")
print(names(sample_read_counts[sample_read_counts < 1e6]))

# Filter out low-depth samples
keep_samples <- sample_read_counts >= 1e6
x <- x[, keep_samples, keep.lib.sizes = FALSE]
dge <- x

###################################################################################################################
# STAGE PROPORTION ANALYSIS
###################################################################################################################
setwd("~/Downloads/")
# Load Su et al. RPKM data and calculate stage proportions
if("Length" %in% colnames(dge$genes) && !"length" %in% colnames(dge$genes)) {
  dge$genes$length <- dge$genes$Length
}
our_log_rpkm <- log2(1 + rpkm(dge))
su_rpkm <- read.csv("./_rpkm_7SexualAndAsexualLifeStages_suetal.csv", header=TRUE, sep=",")
su_rpkm$Ookinete <- NULL; rownames(su_rpkm) <- su_rpkm$ID; su_rpkm$ID <- NULL
su_log_rpkm <- log2(1 + su_rpkm)

# Mixture model fitting function
findMix <- function(Y, X){  
  X[is.na(X)] <- t(replicate(ncol(X), apply(X,1,mean, na.rm=T)))[is.na(X)]
  Rinv <- solve(chol(t(X) %*% X))
  C <- cbind(rep(1,ncol(X)), diag(ncol(X)))
  b <- c(1,rep(0,ncol(X)))
  d <- t(Y) %*% X  
  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq=1)
  sol <- qp$solution; sol[sol<1e-10] <- 0; return(sol)
}

# Fit stage proportions
inter <- intersect(rownames(our_log_rpkm), rownames(su_log_rpkm))
O <- our_log_rpkm[inter, ]; S <- su_log_rpkm[inter, ]
O <- O[order(rownames(O)), ]; S <- S[order(rownames(S)), ]
O <- as.matrix(O); S <- as.matrix(S)
valid_genes <- rownames(O)[complete.cases(O) & complete.cases(S)]
O <- O[valid_genes, ]; S <- S[valid_genes, ]

ourPlotData <- data.frame()
for (i in 1:ncol(O)){
  mix <- findMix(O[,i], as.matrix(S))
  ourPlotData <- rbind(ourPlotData, data.frame(sample=rep(colnames(O)[i], ncol(S)), 
                                               stage=colnames(S), proportion=mix))
}
ourPlotData$stage <- gsub("Gametocyte.*","Gametocyte",ourPlotData$stage)
ourPlotData <- aggregate(proportion~sample+stage,data=ourPlotData,FUN=sum)

# Load sample classifications
classifications_file <- read_excel("classifications_samples.xlsx")
classifications_file <- classifications_file[classifications_file$Sample %in% rownames(dge$samples), ]
dge$samples <- merge(dge$samples, classifications_file, by.x = "row.names", by.y = "Sample", all.x = TRUE)
rownames(dge$samples) <- dge$samples$Row.names; dge$samples$Row.names <- NULL

# Create phenotype categories and add to ourPlotData
categories <- ifelse(dge$samples$Severity == "S", "Severe", "Asymptomatic")
categories <- factor(categories, levels = c("Asymptomatic", "Severe"))
names(categories) <- rownames(dge$samples)
ourPlotData$phenotype <- categories[as.character(ourPlotData$sample)]

###pair-wise comparisons of var domains###
library(dplyr)
library(tidyr)

# 1. Get expression values into a tidy data frame
expr_df <- our_log_rpkm %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID") %>%
  pivot_longer(-GeneID, names_to = "SampleID", values_to = "RPKM")

# 2. Add severity category
categories <- ifelse(dge$samples$Severity == "S", "Severe", "Asymptomatic")
categories <- factor(categories, levels = c("Asymptomatic", "Severe"))
names(categories) <- rownames(dge$samples)

expr_df$Category <- categories[expr_df$SampleID]

# 3. Identify var genes (anything not starting with PF3D7)
expr_df <- expr_df %>%
  mutate(is_var = !grepl("^PF3D7", GeneID))

# -----------------------------
# Total var expression per sample
# -----------------------------
total_var <- expr_df %>%
  filter(is_var) %>%
  group_by(SampleID, Category) %>%
  summarise(TotalVarRPKM = sum(RPKM, na.rm = TRUE), .groups = "drop")

wilcox_total <- wilcox.test(TotalVarRPKM ~ Category, data = total_var)
print(wilcox_total)

# -----------------------------
# Per-domain comparison
# -----------------------------
per_domain <- expr_df %>%
  filter(is_var) %>%
  group_by(SampleID, Category, GeneID) %>%
  summarise(DomainExpr = sum(RPKM, na.rm = TRUE), .groups = "drop")

domain_stats <- per_domain %>%
  group_by(GeneID) %>%
  summarise(
    p_val = wilcox.test(DomainExpr ~ Category)$p.value,
    Severe_mean = mean(DomainExpr[Category == "Severe"]),
    Uncomp_mean = mean(DomainExpr[Category == "Asymptomatic"]),
    .groups = "drop"
  )

domain_stats$padj <- p.adjust(domain_stats$p_val, method = "BH")

write.csv(domain_stats, "Var_Domain_Wilcoxon.csv", row.names = FALSE)

###################################################################################################################
# STAGE PROPORTION PLOT
###################################################################################################################

gg_stages <- ggplot(ourPlotData, aes(x=factor(sample), y=proportion, fill=factor(phenotype))) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("Asymptomatic" = "#d7191c", "Severe" = "#2c7bb6")) +
  facet_wrap(~ stage, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 90),
        axis.text.y = element_text(size=12, angle = 0),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 16, face="bold"),
        legend.text=element_text(size=14),
        legend.key.size = unit(0.25, "in"),
        legend.title = element_text(size=16, face="bold")) +
  labs(x='Sample', y='Proportion') +
  guides(fill=guide_legend(title="Phenotype"))

print(gg_stages)
ggsave("stage_proportions.pdf", gg_stages, width = 16, height = 12)

# Normalization and voom
dge <- calcNormFactors(dge, method="TMM")
dge$samples$group <- categories
design <- model.matrix(~group, data=dge$samples)
v <- voom(dge, design=design, plot=TRUE)

# Stage data preparation
setDT(ourPlotData)
c <- dcast(ourPlotData, sample ~ stage, value.var="proportion")
covs <- data.frame(v$design[,2])
covs <- merge(covs, c, by.x=0, by.y="sample")
colnames(covs) <- c("sample", "disease", colnames(c)[2:ncol(c)])
rownames(covs) <- covs$sample; covs$sample <- NULL
covs <- covs[match(colnames(v$E), rownames(covs)),]

# Control genes
control_genes_gerry <- read.table("./control_genes_gerry.txt", quote="\"")
ctrl_gerry <- geneID_mappings$current[geneID_mappings$current %in% control_genes_gerry$V1]
empirical_controls <- rownames(v$E) %in% ctrl_gerry
if(sum(empirical_controls) < 10) {
  gene_vars <- apply(v$E, 1, var)
  empirical_controls_genes <- names(sort(gene_vars))[1:100]
} else {
  empirical_controls_genes <- rownames(v$E)[empirical_controls]
}

# Analysis setup
FC_THRESHOLD <- 1.0; PVAL_THRESHOLD <- 0.1
output_dir <- "analysis_output"; plots_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results"); gsea_dir <- file.path(output_dir, "gsea_files")
for(dir in c(output_dir, plots_dir, results_dir, gsea_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Disease classifications
severity <- factor(ifelse(grepl("^CH|^ch", colnames(dge)), "Severe", "Asymptomatic"),
                   levels = c("Asymptomatic", "Severe"))
names(severity) <- colnames(dge)
disease_col <- dge$samples$Disease

create_syndrome <- function(pattern) {
  has_syndrome <- grepl(pattern, disease_col) & disease_col != "U"
  is_Asymptomatic <- disease_col == "U"
  result <- ifelse(has_syndrome, TRUE, ifelse(is_Asymptomatic, FALSE, NA))
  names(result) <- colnames(dge); return(result)
}

disease_classifications <- list(
  severity = severity,
  cerebral = create_syndrome("C"),
  hyperparasitemia = create_syndrome("P"), 
  anaemia = create_syndrome("A"),
  hypoglycaemia = create_syndrome("H"),
  acidosis = create_syndrome("L")
)

# Helper functions
compute_real_ruv4 <- function(expression_data, group_factor, control_genes, k_factors = 2) {
  
  # Convert group factor to numeric matrix (like original pipeline)
  categoriesRUV <- data.matrix(as.numeric(group_factor == "Disease"))
  if(length(unique(group_factor)) == 2 && "Severe" %in% levels(group_factor)) {
    categoriesRUV <- data.matrix(as.numeric(group_factor == "Severe"))
  }
  
  # Create empirical controls logical vector
  empirical_controls <- rownames(expression_data) %in% control_genes
  
  # Transpose expression data (genes become rows x samples columns -> samples x genes)
  genes <- data.matrix(t(expression_data))
  
  # Run RUV4 (no Z matrix for now - add staging separately)
  ruv_result <- RUV4(genes, categoriesRUV, empirical_controls, k_factors)
  
  return(list(W = ruv_result$W))
}

create_volcano_plot <- function(limma_results, title, fc_threshold = 1, pval_threshold = 0.1) {
  dge_data <- limma_results %>%
    mutate(gene_type = ifelse(grepl("^PF3D7", gene), "PF3D7", "Other"),
           pval_to_use = pmax(adj.P.Val, .Machine$double.xmin),
           neg_log10_pval = pmin(-log10(pval_to_use), 100),
           significance_category = case_when(
             is.na(adj.P.Val) ~ "NA",
             adj.P.Val < pval_threshold & gene_type == "PF3D7" ~ "PF3D7_Significant",
             adj.P.Val < pval_threshold & gene_type == "Other" ~ "Other_Significant", 
             gene_type == "PF3D7" ~ "PF3D7_Not_Significant",
             TRUE ~ "Other_Not_Significant"))
  
  top_other <- dge_data %>% filter(significance_category == "Other_Significant") %>%
    arrange(adj.P.Val) %>% slice_head(n = 20)
  
  volcano_colors <- c("PF3D7_Significant" = "#FF6347", "Other_Significant" = "#008080", 
                      "PF3D7_Not_Significant" = "#FFB6C1", "Other_Not_Significant" = "#20B2AA", "NA" = "black")
  
  n_pf3d7_sig <- sum(dge_data$significance_category == "PF3D7_Significant", na.rm = TRUE)
  n_other_sig <- sum(dge_data$significance_category == "Other_Significant", na.rm = TRUE)
  y_max <- min(max(dge_data$neg_log10_pval, na.rm = TRUE), 100)
  x_range <- range(dge_data$logFC, na.rm = TRUE); x_max <- max(abs(x_range)) + 1
  
  ggplot(dge_data, aes(x = logFC, y = neg_log10_pval, color = significance_category)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_text_repel(data = top_other, aes(label = gene), size = 2, max.overlaps = 50) +
    scale_color_manual(values = volcano_colors) + theme_minimal() +
    labs(title = title, subtitle = sprintf("PF3D7: %d, Other: %d", n_pf3d7_sig, n_other_sig),
         x = "Log2 Fold Change (Disease vs Asymptomatic)", y = "-Log10 Adjusted P-value") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    xlim(-x_max, x_max) + ylim(0, y_max + 5) + theme(legend.position = "none")
}

# Define colours
colors <- c("Asymptomatic" = "#d7191c", "Severe" = "#2c7bb6")

create_pca_plot <- function(expression_matrix, categories, title) {
  autoplot(
    prcomp(t(expression_matrix)), 
    data = data.frame(Group = categories), 
    colour = 'Group', size = 1, alpha = 0.5
  ) +
    scale_color_manual(values = colors) +
    xlim(-0.15, 0.15) + ylim(-0.15, 0.15) +
    ggtitle(title) + 
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
}


###################################################################################################################
# MA PLOT FUNCTION
###################################################################################################################

library(ggplot2)
library(ggrepel)  # for nice non-overlapping labels

create_ma_plot <- function(limma_results, title) {
  ma_data <- limma_results
  
  # classify genes
  ma_data$gene_type <- ifelse(grepl("^PF3D7", ma_data$gene), "Core", "Var")
  
  # mark significance for var genes only
  ma_data$sig_status <- ifelse(ma_data$gene_type == "Var" & ma_data$adj.P.Val < 0.1,
                               "Var_sig", 
                               ifelse(ma_data$gene_type == "Var", "Var_nonsig", "Core"))
  
  # define colors: Core = grey, Var_sig = dark color, Var_nonsig = light color
  plot_colors <- c("Core" = "grey70",
                   "Var_nonsig" = "#87CEFA",   # light blue
                   "Var_sig" = "#00008B")      # dark blue
  
  ggplot(ma_data, aes(x = AveExpr, y = logFC)) +
    geom_point(aes(color = sig_status), alpha = 0.7, size = 0.7) +
    scale_color_manual(values = plot_colors) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dotted", size = 0.5) +
    # Annotate significant var genes
    geom_text_repel(
      data = subset(ma_data, sig_status == "Var_sig"),
      aes(label = gene),
      size = 2,                 # smaller font
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.5,
      segment.size = 0           # removes lines connecting labels to points
    ) +
    labs(x = "A (Average log2 expression)",
         y = "M (log2 fold change: Disease/Control)",
         title = title) +
    coord_cartesian(xlim = c(-3, 20), ylim = c(-5, 5)) +
    theme_bw() +
    theme(legend.position = "none")
}

get_corrected_expression <- function(original_expr, ruv_factors = NULL, stage_covs = NULL) {
  corrected_expr <- original_expr
  if (!is.null(stage_covs)) {
    for (stage_name in names(stage_covs)) {
      if (is.numeric(stage_covs[[stage_name]])) {
        for (i in 1:nrow(corrected_expr)) {
          lm_fit <- lm(corrected_expr[i, ] ~ stage_covs[[stage_name]])
          corrected_expr[i, ] <- residuals(lm_fit) + mean(corrected_expr[i, ])
        }
      }
    }
  }
  if (!is.null(ruv_factors)) {
    for (i in 1:nrow(corrected_expr)) {
      lm_fit <- lm(corrected_expr[i, ] ~ ruv_factors)
      corrected_expr[i, ] <- residuals(lm_fit) + mean(corrected_expr[i, ])
    }
  }
  return(corrected_expr)
}

create_gct_file <- function(expression_matrix, filename) {
  n_genes <- nrow(expression_matrix); n_samples <- ncol(expression_matrix)
  gct_lines <- character()
  gct_lines[1] <- "#1.2"; gct_lines[2] <- paste(n_genes, n_samples, sep = "\t")
  header <- paste(c("NAME", "DESCRIPTION", colnames(expression_matrix)), collapse = "\t")
  gct_lines[3] <- header
  gene_names <- rownames(expression_matrix)
  for (i in 1:n_genes) {
    data_line <- paste(c(gene_names[i], gene_names[i], expression_matrix[i, ]), collapse = "\t")
    gct_lines[3 + i] <- data_line
  }
  writeLines(gct_lines, filename)
}

create_cls_file <- function(sample_groups, filename, class_names = NULL) {
  n_samples <- length(sample_groups)
  unique_groups <- unique(sample_groups); n_classes <- length(unique_groups)
  if (is.null(class_names)) class_names <- unique_groups
  numeric_assignments <- match(sample_groups, unique_groups) - 1
  cls_lines <- character()
  cls_lines[1] <- paste(n_samples, n_classes, "1", sep = " ")
  cls_lines[2] <- paste("#", paste(class_names, collapse = " "))
  cls_lines[3] <- paste(numeric_assignments, collapse = " ")
  writeLines(cls_lines, filename)
}

create_gsea_files <- function(expression_matrix, sample_groups, method_name, disease_name, output_dir) {
  clean_method <- gsub("[^A-Za-z0-9_]", "_", method_name)
  clean_disease <- gsub("[^A-Za-z0-9_]", "_", disease_name)
  base_filename <- paste(clean_method, clean_disease, sep = "_")
  gct_filename <- file.path(output_dir, paste0(base_filename, ".gct"))
  cls_filename <- file.path(output_dir, paste0(base_filename, ".cls"))
  create_gct_file(expression_matrix, gct_filename)
  create_cls_file(as.character(sample_groups), cls_filename)
  return(list(gct = gct_filename, cls = cls_filename))
}

# Stage data preparation
setDT(ourPlotData)
c <- dcast(ourPlotData, sample ~ stage, value.var = "proportion")
covs_base <- data.frame(row.names = colnames(v$E))
covs_base <- merge(covs_base, c, by.x = 0, by.y = "sample")
colnames(covs_base) <- c("sample", colnames(c)[2:ncol(c)])
rownames(covs_base) <- covs_base$sample; covs_base$sample <- NULL
covs_base <- covs_base[match(colnames(v$E), rownames(covs_base)), ]

all_results <- list(); all_plots <- list()

# Main analysis loop
for(disease_name in names(disease_classifications)) {
  selected_class <- disease_classifications[[disease_name]]
  selected_class_vector <- as.vector(unlist(selected_class))
  
  if(is.logical(selected_class_vector)) {
    valid <- !is.na(selected_class_vector)
    filtered_expression <- v$E[, valid, drop=FALSE]
    filtered_covs <- covs_base[valid, , drop=FALSE]
    grp <- factor(ifelse(selected_class_vector[valid], "Disease", "Asymptomatic"),
                  levels = c("Asymptomatic", "Disease"))
  } else if(disease_name == "severity") {
    valid <- !is.na(selected_class_vector)
    filtered_expression <- v$E[, valid, drop=FALSE]
    filtered_covs <- covs_base[valid, , drop=FALSE]
    grp <- selected_class[valid]
  }
  
  covs_use <- filtered_covs
  for(stage_col in c("Ring", "Gametocyte", "Schizont")) {
    if(stage_col %in% colnames(covs_use) && is.numeric(covs_use[[stage_col]])) {
      covs_use[[stage_col]] <- scale(covs_use[[stage_col]], center = TRUE, scale = FALSE)
    }
  }
  
  df_base <- data.frame(disease = grp, covs_use)
  D0 <- model.matrix(~ 0 + disease, data = df_base)
  D_ring <- cbind(D0, Ring = df_base$Ring)
  D_ring_all <- cbind(D0, Ring = df_base$Ring, Gametocyte = df_base$Gametocyte, Schizont = df_base$Schizont)
  
  colnames(D0) <- make.names(colnames(D0))
  colnames(D_ring) <- make.names(colnames(D_ring))
  colnames(D_ring_all) <- make.names(colnames(D_ring_all))
  
  ruv_factors <- compute_real_ruv4(filtered_expression, grp, empirical_controls_genes, k_factors = 2)
  D_ruv_no <- cbind(D0, ruv_factors$W)
  D_ring_ruv <- cbind(D_ring, ruv_factors$W)
  D_all_ruv <- cbind(D_ring_all, ruv_factors$W)
  
  design_list <- list(D0, D_ruv_no, D_ring, D_ring_ruv, D_ring_all, D_all_ruv)
  method_names <- c("Library_Size_Only", "Library_Size_RUV", "Library_Size_Ring", 
                    "Library_Size_Ring_RUV", "Library_Size_All_Stages", "Library_Size_All_Stages_RUV")
  
  limma_results_list <- list(); volcano_plots <- list(); pca_plots <- list(); ma_plots <- list()
  
  for(i in seq_along(design_list)) {
    design <- design_list[[i]]
    fit <- limma::lmFit(filtered_expression, design)
    
    design_cols <- colnames(design)
    if(disease_name == "severity") {
      disease_col_fixed <- design_cols[grepl("Severe", design_cols)][1]
      control_col_fixed <- design_cols[grepl("Asymptomatic", design_cols)][1]
    } else {
      disease_col_fixed <- design_cols[grepl("Disease", design_cols)][1]
      control_col_fixed <- design_cols[grepl("Asymptomatic", design_cols)][1]
    }
    
    if(!is.na(disease_col_fixed) && !is.na(control_col_fixed)) {
      contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
      rownames(contrast_matrix) <- colnames(design)
      colnames(contrast_matrix) <- "Disease_vs_Control"
      contrast_matrix[disease_col_fixed, 1] <- 1
      contrast_matrix[control_col_fixed, 1] <- -1
      fit2 <- limma::contrasts.fit(fit, contrast_matrix)
      fit2 <- limma::eBayes(fit2)
      results <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    } else {
      fit2 <- limma::eBayes(fit)
      coef_idx <- which(grepl("Disease|Severe", colnames(fit2$coefficients)))[1]
      results <- limma::topTable(fit2, coef = ifelse(!is.na(coef_idx), coef_idx, 2), 
                                 number = Inf, sort.by = "none")
    }
    
    results$gene <- rownames(results)
    limma_results_list[[i]] <- results
    
    # Create plots
    volcano_title <- paste("Volcano:", gsub("_", " ", method_names[i]), "-", disease_name)
    volcano_plots[[i]] <- create_volcano_plot(results, volcano_title, FC_THRESHOLD, PVAL_THRESHOLD)
    
    # Corrected expression for visualization
    stage_covs_list <- NULL; ruv_use <- NULL
    if (grepl("Ring", method_names[i]) && !grepl("All_Stages", method_names[i])) {
      stage_covs_list <- list(Ring = df_base$Ring)
    } else if (grepl("All_Stages", method_names[i])) {
      stage_covs_list <- list(Ring = df_base$Ring, Gametocyte = df_base$Gametocyte, Schizont = df_base$Schizont)
    }
    if (grepl("RUV", method_names[i])) ruv_use <- ruv_factors$W
    
    corrected_expr <- get_corrected_expression(filtered_expression, ruv_use, stage_covs_list)
    pca_plots[[i]] <- create_pca_plot(corrected_expr, grp, paste("PCA:", gsub("_", " ", method_names[i])))
    ma_plots[[i]] <- create_ma_plot(results, paste("MA:", gsub("_", " ", method_names[i])))
    
    # Save individual MA plot
    ma_plot_filename <- file.path(plots_dir, paste0("MA_", method_names[i], "_", disease_name, ".pdf"))
    ggsave(filename = ma_plot_filename, plot = ma_plots[[i]], width = 10, height = 8)
    
    # Save results and GSEA files
    write.csv(results, file.path(results_dir, paste0(method_names[i], "_", disease_name, "_results.csv")))
    create_gsea_files(corrected_expr, grp, method_names[i], disease_name, gsea_dir)
  }
  
  # Save grid plots
  pdf(file.path(plots_dir, paste0("PCA_grid_", disease_name, ".pdf")), width = 12, height = 12)
  grid.arrange(grobs = pca_plots, ncol = 2, nrow = 3, top = paste("PCA Comparison:", disease_name))
  dev.off()
  
  pdf(file.path(plots_dir, paste0("MA_grid_", disease_name, ".pdf")), width = 30, height = 30)
  grid.arrange(grobs = ma_plots, ncol = 2, nrow = 3, top = paste("MA Comparison:", disease_name))
  dev.off()
  
  pdf(file.path(plots_dir, paste0("Volcano_grid_", disease_name, ".pdf")), width = 12, height = 12)
  grid.arrange(grobs = volcano_plots, ncol = 2, nrow = 3, top = paste("Volcano Comparison:", disease_name))
  dev.off()
  
  all_results[[disease_name]] <- limma_results_list
  all_plots[[disease_name]] <- volcano_plots
}

# =============================================================================
# GENE ANNOTATION FUNCTIONS
# =============================================================================

# Function to read gene annotations from Excel file
read_gene_annotations <- function(file_path) {
  if (!file.exists(file_path)) {
    cat("Warning: Gene annotations file not found:", file_path, "\n")
    return(NULL)
  }
  
  tryCatch({
    annotations <- read_excel(file_path)
    cat("Successfully loaded", nrow(annotations), "gene annotations\n")
    return(annotations)
  }, error = function(e) {
    cat("Error reading gene annotations:", e$message, "\n")
    return(NULL)
  })
}

# Function to export DE results with annotations
export_de_results <- function(limma_results, method_name, analysis_type, 
                              fc_threshold = 1.0, pval_threshold = 0.05,
                              gene_annotations = NULL) {
  
  # Add gene column if not present
  if (!"gene" %in% colnames(limma_results)) {
    limma_results$gene <- rownames(limma_results)
  }
  
  # Add annotations if provided
  if (!is.null(gene_annotations)) {
    # Assume the gene annotations have a column that matches with gene IDs
    # You may need to adjust the column names based on your Excel file structure
    annotation_cols <- colnames(gene_annotations)
    gene_id_col <- annotation_cols[1]  # Assume first column contains gene IDs
    
    limma_results <- merge(limma_results, gene_annotations, 
                          by.x = "gene", by.y = gene_id_col, 
                          all.x = TRUE, sort = FALSE)
  }
  
  # Add significance flags
  limma_results$significant <- !is.na(limma_results$adj.P.Val) & 
                              limma_results$adj.P.Val < pval_threshold
  
  limma_results$fc_threshold <- !is.na(limma_results$logFC) & 
                               abs(limma_results$logFC) >= fc_threshold
  
  limma_results$both_significant <- limma_results$significant & limma_results$fc_threshold
  
  # Sort by adjusted p-value
  limma_results <- limma_results[order(limma_results$adj.P.Val, na.last = TRUE), ]
  
  # Create filename
  filename <- paste0(analysis_type, "_", method_name, "_annotated_results.csv")
  filepath <- file.path(results_dir, filename)
  
  # Export results
  write.csv(limma_results, filepath, row.names = FALSE)
  cat("Exported", method_name, "results to:", filename, "\n")
  cat("  - Total genes:", nrow(limma_results), "\n")
  cat("  - Significant (adj.P.Val <", pval_threshold, "):", sum(limma_results$significant, na.rm = TRUE), "\n")
  cat("  - FC threshold (|logFC| >=", fc_threshold, "):", sum(limma_results$fc_threshold, na.rm = TRUE), "\n")
  cat("  - Both significant:", sum(limma_results$both_significant, na.rm = TRUE), "\n")
  
  return(limma_results)
}

# =============================================================================
# LOAD GENE ANNOTATIONS
# =============================================================================
cat("=== LOADING GENE ANNOTATIONS ===\n")
gene_annotations <- read_gene_annotations("gene_names.xlsx")

# =============================================================================
# EXPORT ALL LIMMA-VOOM RESULTS WITH ANNOTATIONS
# =============================================================================
cat("=== EXPORTING LIMMA-VOOM RESULTS WITH ANNOTATIONS ===\n")

# Export all results with annotations for each disease classification
limma_exported <- list()

for(disease_name in names(all_results)) {
  cat("Processing disease classification:", disease_name, "\n")
  
  # Method names from your existing analysis
  method_names <- c("Library_Size_Only", "Library_Size_RUV", "Library_Size_Ring", 
                    "Library_Size_Ring_RUV", "Library_Size_All_Stages", "Library_Size_All_Stages_RUV")
  
  limma_exported[[disease_name]] <- list()
  
  for(i in seq_along(method_names)) {
    method_name <- method_names[i]
    limma_results <- all_results[[disease_name]][[i]]
    
    if (!is.null(limma_results)) {
      limma_exported[[disease_name]][[method_name]] <- export_de_results(
        limma_results, 
        paste(method_name, disease_name, sep = "_"), 
        "limma", 
        fc_threshold = FC_THRESHOLD, 
        pval_threshold = PVAL_THRESHOLD,
        gene_annotations = gene_annotations
      )
    }
  }
}

# =============================================================================
# CREATE SUMMARY TABLES
# =============================================================================
cat("=== CREATING SUMMARY TABLES ===\n")

# Create summary of significant genes across methods and diseases
summary_results <- data.frame()

for(disease_name in names(limma_exported)) {
  for(method_name in names(limma_exported[[disease_name]])) {
    results <- limma_exported[[disease_name]][[method_name]]
    if (!is.null(results)) {
      summary_row <- data.frame(
        Disease = disease_name,
        Method = method_name,
        Total_Genes = nrow(results),
        Significant_Genes = sum(results$significant, na.rm = TRUE),
        FC_Threshold_Genes = sum(results$fc_threshold, na.rm = TRUE),
        Both_Significant = sum(results$both_significant, na.rm = TRUE),
        Core_Genes_Significant = sum(results$significant & grepl("^PF3D7", results$gene), na.rm = TRUE),
        Var_Genes_Significant = sum(results$significant & !grepl("^PF3D7", results$gene), na.rm = TRUE)
      )
      summary_results <- rbind(summary_results, summary_row)
    }
  }
}

# Export summary
write.csv(summary_results, file.path(results_dir, "limma_analysis_summary.csv"), row.names = FALSE)
cat("Summary table exported to: limma_analysis_summary.csv\n")

# Print summary
print(summary_results)
###################################################################################################################
# COMPREHENSIVE HEATMAP ANALYSIS WITH CORE AND VAR GENE SEPARATION
###################################################################################################################

# Heatmap for severity analysis
if("severity" %in% names(all_results)) {  
  method_index <- 6  # Library_Size_All_Stages_RUV
  limma_results <- all_results[["severity"]][[method_index]]
  significant_genes <- limma_results[!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1, ]
  
  if(nrow(significant_genes) > 0) {    
    selected_class <- disease_classifications[["severity"]]
    valid <- !is.na(selected_class)
    filtered_expression <- v$E[, valid, drop=FALSE]
    
    gene_ids <- rownames(significant_genes)
    matching_genes <- intersect(rownames(filtered_expression), gene_ids)
    heatmap_data <- filtered_expression[matching_genes, , drop = FALSE]
    
    # Separate core and var genes
    core_genes <- matching_genes[grepl("^PF3D7", matching_genes)]
    var_genes <- matching_genes[!grepl("^PF3D7", matching_genes)]
    
    # Sample annotations
    disease_codes <- dge$samples[colnames(heatmap_data), "Disease"]
    disease_labels <- case_when(
      disease_codes == "U" ~ "Asymptomatic",
      str_detect(disease_codes, "C") ~ "Cerebral",
      (str_detect(disease_codes, "L") | str_detect(disease_codes, "H")) & !str_detect(disease_codes, "C") ~ "Acidosis_Hypo",
      TRUE ~ "Other"
    )
    
    col_annot <- data.frame(Phenotype = factor(disease_labels), row.names = colnames(heatmap_data))
    disease_colors <- c("Asymptomatic" = "#d7191c", "Cerebral" = "slateblue4", 
                        "Acidosis_Hypo" = "deepskyblue2", "Other" = "turquoise1")
    annotation_colors <- list(Phenotype = disease_colors)
    
    # Create separate heatmaps for core and var genes
    if(length(core_genes) > 0) {  
      # Save the pheatmap object
      ph <- pheatmap(
        heatmap_data[core_genes, , drop=FALSE], 
        show_rownames = ifelse(length(core_genes) <= 60, TRUE, FALSE),
        show_colnames = FALSE, scale = "row", 
        clustering_distance_rows = "correlation",
        cluster_cols = TRUE,
        annotation_col = col_annot, 
        annotation_colors = annotation_colors,
        breaks = seq(-2, 2, length.out = 100),
        main = paste("DE Core Genes: Library Size + All Stages + RUV\nSeverity classification"),
        filename = file.path(plots_dir, "heatmap_severity_core_genes.pdf"),
        width = 20, height = min(16, max(10, length(core_genes) * 0.15 + 6))
      )
     
      # After creating your pheatmap object 'ph':
      sample_tree <- ph$tree_col
      
      # 1. Get original dendrogram clusters (11 clusters)
      numeric_clusters <- cutree(sample_tree, k = 11)
      
      # 2. Identify asymptomatic and severe samples
      all_samples <- colnames(filtered_expression)
      
      # disease_labels is a positional vector matching filtered_expression columns
      asymp_indices <- which(disease_labels == "Asymptomatic")
      severe_indices <- which(disease_labels != "Asymptomatic")
      
      asymp_samples <- all_samples[asymp_indices]
      severe_samples <- all_samples[severe_indices]
      
      cat("Found", length(asymp_samples), "asymptomatic samples\n")
      cat("Found", length(severe_samples), "severe samples\n")
      
      # 3. Create new cluster assignments with simple names
      cluster_labels <- rep(NA, length(all_samples))
      names(cluster_labels) <- all_samples
      
      # Assign severe samples to simple numbered clusters (S1, S2, etc.)
      for(sample in severe_samples) {
        if(sample %in% names(numeric_clusters)) {
          cluster_labels[sample] <- paste0("S", numeric_clusters[sample])
        }
      }
      
      # Assign ALL asymptomatic samples to "ASYMP"
      cluster_labels[asymp_samples] <- "ASYMP"
      
      # 4. Remove any samples that weren't assigned (shouldn't happen but safety check)
      valid_samples <- !is.na(cluster_labels)
      cluster_labels <- cluster_labels[valid_samples]
      gct_matrix <- as.matrix(filtered_expression[, names(cluster_labels)])
      
      # 5. Verify we have the right structure
      cat("Cluster summary:\n")
      print(table(cluster_labels))
      cat("\nNumber of unique clusters:", length(unique(cluster_labels)), "\n")
      
      # 6. Prepare GCT file
      num_genes <- nrow(gct_matrix)
      num_samples <- ncol(gct_matrix)
      gct_header <- c("#1.2", paste(num_genes, num_samples, sep="\t"))
      
      # Add gene descriptions (use gene names as descriptions)
      gene_descriptions <- rownames(gct_matrix)
      gct_content <- cbind(rownames(gct_matrix), gene_descriptions, gct_matrix)
      
      gct_file <- c(
        gct_header,
        paste(c("Name", "Description", colnames(gct_matrix)), collapse="\t"),
        apply(gct_content, 1, function(x) paste(x, collapse="\t"))
      )
      
      # 7. Prepare CLS file
      unique_classes <- unique(cluster_labels)
      num_classes <- length(unique_classes)
      
      # Debug: check cluster_labels before creating CLS
      cat("Cluster labels summary before CLS creation:\n")
      print(table(cluster_labels))
      cat("First 10 cluster labels:", head(cluster_labels, 10), "\n")
      cat("Last 10 cluster labels:", tail(cluster_labels, 10), "\n")
      
      cls_file <- c(
        paste(num_samples, num_classes, 1),
        paste("#", paste(unique_classes, collapse=" ")),
        paste(cluster_labels, collapse=" ")
      )
      
      # 8. Save to Desktop
      desktop_path <- "~/Desktop"
      writeLines(gct_file, file.path(desktop_path, "severity_clusters.gct"))
      writeLines(cls_file, file.path(desktop_path, "severity_clusters.cls"))
      
      # 9. Create a summary file showing which samples are in which clusters
      cluster_summary <- data.frame(
        Sample = names(cluster_labels),
        Cluster = cluster_labels,
        Original_Phenotype = disease_labels[match(names(cluster_labels), all_samples)],
        stringsAsFactors = FALSE
      )
      write.csv(cluster_summary, file.path(desktop_path, "cluster_assignments.csv"), row.names = FALSE)
      
      message("GCT + CLS files created successfully!")
      message("- ", sum(cluster_labels == "ASYMP"), " asymptomatic samples in cluster 'ASYMP'")
      message("- ", sum(cluster_labels != "ASYMP"), " severe samples distributed across ", 
              length(unique(cluster_labels[cluster_labels != "ASYMP"])), " severe clusters (S1, S2, S3, etc.)")
    }    
  }
}
    #sig var genes only#
    if(length(var_genes) > 0) {
      pheatmap(heatmap_data[var_genes, , drop=FALSE], 
               show_rownames = ifelse(length(var_genes) <= 60, TRUE, FALSE),
               show_colnames = FALSE, scale = "row", clustering_distance_rows = "correlation",
               cluster_cols = TRUE, annotation_col = col_annot, annotation_colors = annotation_colors,
               breaks = seq(-2, 2, length.out = 100),
               main = paste("DE Var Genes: Library Size + All Stages + RUV\nSeverity classification"),
               filename = file.path(plots_dir, "heatmap_severity_sig_var_genes.pdf"),
               width = 20, height = min(16, max(10, length(var_genes) * 0.15 + 6)))
    }
    
    # Combined heatmap
    pheatmap(heatmap_data, show_rownames = ifelse(length(matching_genes) <= 60, TRUE, FALSE),
             show_colnames = FALSE, scale = "row", clustering_distance_rows = "correlation",
             cluster_cols = TRUE, annotation_col = col_annot, annotation_colors = annotation_colors,
             breaks = seq(-2, 2, length.out = 100),
             main = paste("DE Genes: Library Size + All Stages + RUV\nSeverity classification"),
             filename = file.path(plots_dir, "heatmap_severity_all_genes.pdf"),
             width = 20, height = min(16, max(10, length(matching_genes) * 0.15 + 6)))

###for all var genes###

if("severity" %in% names(all_results)) {
  method_index <- 6  # Library_Size_All_Stages_RUV
  limma_results <- all_results[["severity"]][[method_index]]
  significant_genes <- limma_results[!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1, ]
  
  if(nrow(significant_genes) > 0) {
    selected_class <- disease_classifications[["severity"]]
    valid <- !is.na(selected_class)
    filtered_expression <- v$E[, valid, drop=FALSE]
    
    # For core genes - use significant genes only
    gene_ids <- rownames(significant_genes)
    matching_genes <- intersect(rownames(filtered_expression), gene_ids)
    heatmap_data <- filtered_expression[matching_genes, , drop = FALSE]
    
    # Separate core and var genes
    core_genes <- matching_genes[grepl("^PF3D7", matching_genes)]
    var_genes_sig <- matching_genes[!grepl("^PF3D7", matching_genes)]
    
    # For var genes - use all var genes (no significance filtering)
    all_var_genes <- rownames(filtered_expression)[!grepl("^PF3D7", rownames(filtered_expression))]
    var_heatmap_data <- filtered_expression[all_var_genes, , drop = FALSE]
    
    # Sample annotations
    disease_codes <- dge$samples[colnames(filtered_expression), "Disease"]
    disease_labels <- case_when(
      disease_codes == "U" ~ "Asymptomatic",
      str_detect(disease_codes, "C") ~ "Cerebral",
      (str_detect(disease_codes, "L") | str_detect(disease_codes, "H")) & !str_detect(disease_codes, "C") ~ "Acidosis_Hypo",
      TRUE ~ "Other"
    )
    
    col_annot <- data.frame(Phenotype = factor(disease_labels), row.names = colnames(filtered_expression))
    disease_colors <- c("Asymptomatic" = "#d7191c", "Cerebral" = "slateblue4", 
                        "Acidosis_Hypo" = "deepskyblue2", "Other" = "turquoise1")
    annotation_colors <- list(Phenotype = disease_colors)

    if(length(all_var_genes) > 0) {
      pheatmap(var_heatmap_data, 
               show_rownames = ifelse(length(all_var_genes) <= 160, TRUE, FALSE),
               show_colnames = FALSE, scale = "row", clustering_distance_rows = "correlation",
               cluster_cols = TRUE, annotation_col = col_annot, annotation_colors = annotation_colors,
               breaks = seq(-2, 2, length.out = 100),
               main = paste("All Var Genes: Library Size + All Stages + RUV\nSeverity classification"),
               filename = file.path(plots_dir, "heatmap_severity_all_var_genes.pdf"),
               width = 20, height = min(30, max(20, length(all_var_genes) * 0.15 + 6)))
    }

  }
}

###RPKM based heatmap###

if("severity" %in% names(all_results)) {
  method_index <- 6  # Library_Size_All_Stages_RUV
  limma_results <- all_results[["severity"]][[method_index]]
  
  # Significant genes: adj.P.Val < 0.1
  significant_genes <- limma_results[!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1, ]
  sig_gene_ids <- rownames(significant_genes)
  
  # Filter RPKM data to only var genes
  var_genes_all <- rownames(our_log_rpkm)[!grepl("^PF3D7", rownames(our_log_rpkm))]
  var_genes_sig <- intersect(var_genes_all, sig_gene_ids)
  
  # FIXED: Only use samples that exist in both our_log_rpkm AND dge$samples
  common_samples <- intersect(colnames(our_log_rpkm), rownames(dge$samples))
  our_log_rpkm_filtered <- our_log_rpkm[, common_samples, drop=FALSE]
  
  # Prepare sample annotations (syndromes) - only for common samples
  disease_codes <- dge$samples[common_samples, "Disease"]
  disease_labels <- case_when(
    disease_codes == "U" ~ "Asymptomatic",
    str_detect(disease_codes, "C") ~ "Cerebral",
    (str_detect(disease_codes, "L") | str_detect(disease_codes, "H")) & !str_detect(disease_codes, "C") ~ "Acidosis_Hypo",
    TRUE ~ "Other"
  )
  col_annot <- data.frame(Phenotype = factor(disease_labels), row.names = common_samples)
  disease_colors <- c(
    "Asymptomatic" = "#d7191c", 
    "Cerebral" = "slateblue4", 
    "Acidosis_Hypo" = "deepskyblue2", 
    "Other" = "turquoise1"
  )
  annotation_colors <- list(Phenotype = disease_colors)
  
  # Heatmap for ALL var genes
  heatmap_all_var <- our_log_rpkm_filtered[var_genes_all, , drop=FALSE]
  
  pheatmap(heatmap_all_var,
           show_rownames = TRUE,
           show_colnames = FALSE,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = col_annot,
           annotation_colors = annotation_colors,
           main = "All Var Genes (log2 RPKM) with Syndromes",
           filename = file.path(plots_dir, "heatmap_rpkm_all_var_genes.pdf"),
           width = 20,
           height = min(50, max(45, nrow(heatmap_all_var) * 0.15 + 6)))
  
  # Heatmap for SIGNIFICANT var genes
  if(length(var_genes_sig) > 0) {
    heatmap_sig_var <- our_log_rpkm_filtered[var_genes_sig, , drop=FALSE]
    
    pheatmap(heatmap_sig_var,
             show_rownames = TRUE,
             show_colnames = FALSE,
             scale = "none",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_col = col_annot,
             annotation_colors = annotation_colors,
             main = paste0("Significant Var Genes (adj.P.Val < 0.1, log2 RPKM), n=", length(var_genes_sig)),
             filename = file.path(plots_dir, "heatmap_rpkm_sig_var_genes.pdf"),
             width = 20,
             height = min(16, max(10, nrow(heatmap_sig_var) * 0.15 + 6)))
  } else {
    cat("No significant var genes found (adj.P.Val < 0.1)\n")
  }
}

####bootstrapped dendrogram for sig var genes####

if(length(var_genes) > 1) {
  library(pvclust)
  
  # Run pvclust
  pv_result <- pvclust(
    t(heatmap_data[var_genes, , drop = FALSE]),  # transpose if needed
    method.hclust = "average",
    method.dist = "correlation",
    nboot = 1000,
    parallel = TRUE
  )
  
  # Save dendrogram
  pdf(file.path(plots_dir, "pvclust_dendrogram_var_genes.pdf"), width=20, height=10)
  plot(pv_result, main="Bootstrap dendrogram (var genes)")
  pvrect(pv_result, alpha=0.90)
  dev.off()
  
  # Use order from pvclust
  ordered_genes <- rownames(heatmap_data[var_genes, ])[pv_result$hclust$order]
} else {
  ordered_genes <- var_genes
}

###bootstrapped dendogram for all var genes###

if(length(rownames(heatmap_data)) > 1) {
  library(pvclust)
  
  # ⚡ Redefine var_genes as ALL var genes in v$E (not filtered by significance)
  all_var_genes <- rownames(v$E)[!grepl("^PF3D7", rownames(v$E))]
  all_var_genes <- intersect(all_var_genes, rownames(v$E))  # safety check
  
  # Make expression matrix for all var genes
  heatmap_data_all <- v$E[all_var_genes, , drop = FALSE]
  
  # Run pvclust
  pv_result <- pvclust(
    t(heatmap_data_all),          # transpose: samples × genes
    method.hclust = "average",
    method.dist = "correlation",
    nboot = 1000,
    parallel = TRUE
  )
  
  # Save dendrogram
  pdf(file.path(plots_dir, "pvclust_dendrogram_all_var_genes.pdf"), width=40, height=10)
  plot(pv_result, main="Bootstrap dendrogram (all var genes)")
  pvrect(pv_result, alpha=0.90)
  dev.off()
  
  # Use order from pvclust
  ordered_genes <- rownames(heatmap_data_all)[pv_result$hclust$order]
} else {
  ordered_genes <- var_genes
}

###################################################################################################################
# BOXPLOTS FOR VAR GENES
###################################################################################################################
# Boxplots for var genes
var_genes <- rownames(v$E)[!grepl("^PF3D7", rownames(v$E))]
var_expr <- v$E[var_genes, ]

# Log CPM boxplots
var_df <- as.data.frame(var_expr)
var_df$GeneID <- rownames(var_expr)
var_long <- var_df %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  mutate(severity = severity[sample])

n_genes_per_page <- 53
n_pages <- ceiling(length(var_genes) / n_genes_per_page)
for(page in 1:n_pages){
  genes_this_page <- var_genes[((page-1)*n_genes_per_page + 1) : min(page*n_genes_per_page, length(var_genes))]
  plot_page <- var_long %>%
    filter(GeneID %in% genes_this_page) %>%
    ggplot(aes(x = GeneID, y = expression, fill = severity)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("Asymptomatic" = "#d7191c", "Severe" = "#2c7bb6")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Var gene expression by severity (Page", page, "of", n_pages, ")"),
         x = "Var gene / domain", y = "Expression")
  ggsave(filename = file.path(plots_dir, paste0("var_boxplots_page", page, ".pdf")), 
         plot = plot_page, width = 16, height = 8)
}

# Z-score normalized boxplots
var_expr_z <- t(apply(var_expr, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
var_df_z <- as.data.frame(var_expr_z)
var_df_z$GeneID <- rownames(var_expr_z)
var_long_z <- var_df_z %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  mutate(severity = severity[sample])

for(page in 1:n_pages){
  genes_this_page <- var_genes[((page-1)*n_genes_per_page + 1) : min(page*n_genes_per_page, length(var_genes))]
  plot_page <- var_long_z %>%
    filter(GeneID %in% genes_this_page) %>%
    ggplot(aes(x = GeneID, y = expression, fill = severity)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("Asymptomatic" = "#d7191c", "Severe" = "#2c7bb6")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Var gene expression by severity (Page", page, "of", n_pages, ")"),
         x = "Var gene / domain", y = "Expression (z-score)")
  ggsave(filename = file.path(plots_dir, paste0("var_boxplots_z_page", page, ".pdf")), 
         plot = plot_page, width = 16, height = 8)
}

###################################################################################################################
# ICAM1 binding motif boxplots
###################################################################################################################
####plotting based on cerebral malaria or severe malaria but not cerebral####
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
# Subset to ICAM1-binding DBLb domains
icam1_domains <- c("DBLb1","DBLb11","DBLb12","DBLb3","DBLb4","DBLb5","DBLb8")

icam_hits <- read_excel("final_filtered_icam_binding.xlsx")

# Only keep samples that had ICAM1 hits
icam1_samples <- unique(icam_hits$Sample)

# Intersect with samples that have CM info
valid_samples <- intersect(icam1_samples, names(disease_classifications$cerebral))

# Subset expression matrix BEFORE pivoting
dblb_expr <- v$E[rownames(v$E) %in% icam1_domains, valid_samples]

# Convert to long format
dblb_long <- as.data.frame(dblb_expr) %>%
  tibble::rownames_to_column("GeneID") %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  mutate(CM_status = ifelse(disease_classifications$cerebral[sample], "severe_cerebral", "severe_nonCerebral")) %>%
  filter(!is.na(CM_status))  # <<< REMOVE rows with NA

# Aggregate across domains per sample
dblb_agg <- dblb_long %>%
  group_by(sample, CM_status) %>%
  summarise(expression_sum = sum(expression), .groups = "drop")

# Plot
ggplot(dblb_agg, aes(x = CM_status, y = expression_sum, fill = CM_status)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("severe_cerebral" = "#d7191c", "severe_nonCerebral" = "#2c7bb6")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) +
  labs(title = "ICAM1-binding DBLb domain expression",
       x = "Cerebral Malaria status", y = "Aggregated logCPM expression")

####plotting based on DC4 vs non-DC4 ICAM-1 binding####
# -------------------------
# Parameters
# -------------------------
icam1_domains <- c("DBLb1","DBLb11","DBLb12",
                   "DBLb3","DBLb4","DBLb5","DBLb8")

# Load FIMO-positive ICAM1 hits
icam_hits <- read_excel("final_filtered_icam_binding.xlsx")

# DC4 samples: motif_alt_id == "ICAM1_1"
dc4_samples <- unique(icam_hits$Sample[icam_hits$motif_alt_id == "ICAM1_1"])

# -------------------------
# Subset expression matrix
# -------------------------
dblb_expr <- v$E[rownames(v$E) %in% icam1_domains, ]

# Keep only samples that have CM info (to be consistent with previous scripts)
samples_to_keep <- intersect(colnames(dblb_expr), names(disease_classifications$cerebral))
dblb_expr <- dblb_expr[, samples_to_keep]

# -------------------------
# Prepare long-format dataframe
# -------------------------
dblb_df <- as.data.frame(dblb_expr)
dblb_df$GeneID <- rownames(dblb_expr)

dblb_long <- dblb_df %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  mutate(
    DC4_status = ifelse(sample %in% dc4_samples, "DC4", "Non-DC4")
  ) %>%
  filter(!is.na(DC4_status)) %>%
  mutate(DC4_status = factor(DC4_status, levels = c("Non-DC4", "DC4")))

# -------------------------
# Aggregate across domains per sample
# -------------------------
dblb_agg <- dblb_long %>%
  group_by(sample, DC4_status) %>%
  summarise(expression_sum = sum(expression), .groups = "drop")

# -------------------------
# Plot boxplot
# -------------------------
plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

ggplot(dblb_agg, aes(x = DC4_status, y = expression_sum, fill = DC4_status)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("DC4" = "#d7191c", "Non-DC4" = "#2c7bb6")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) +
  labs(title = "DBLb domain expression by DC4 status",
       x = "DC4 status", y = "Aggregated logCPM expression")

ggsave(file.path(plots_dir, "DBLb_DC4_boxplot.pdf"), width = 6, height = 6)

# -------------------------
# Add CM status back in
# -------------------------
dblb_long <- dblb_df %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  mutate(
    CM_status = ifelse(disease_classifications$cerebral[sample],
                       "severe_cerebral", "severe_nonCerebral"),
    DC4_status = ifelse(sample %in% dc4_samples, "DC4", "Non-DC4")
  ) %>%
  filter(!is.na(CM_status), !is.na(DC4_status))

dblb_agg <- dblb_long %>%
  group_by(sample, CM_status, DC4_status) %>%
  summarise(expression_sum = sum(expression), .groups = "drop")

# -------------------------
# Wilcoxon tests
# -------------------------
# Non-DC4
wilcox_nonDC4 <- wilcox.test(
  expression_sum ~ CM_status,
  data = filter(dblb_agg, DC4_status == "Non-DC4"),
  alternative = "greater"  # tests if CM > non-CM
)

# DC4
wilcox_DC4 <- wilcox.test(
  expression_sum ~ CM_status,
  data = filter(dblb_agg, DC4_status == "DC4"),
  alternative = "greater"  # tests if CM > non-CM
)

wilcox_nonDC4
wilcox_DC4


####plotting both cerebral malaria status and DC4 status####
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# -------------------------
# Parameters
# -------------------------
icam1_domains <- c("DBLb1","DBLb11","DBLb12",
                   "DBLb3","DBLb4","DBLb5","DBLb8")

# Load FIMO-positive ICAM1 hits
# 'final_filtered_icam_binding.xlsx' has columns: Sample, sequence_name, motif_alt_id
icam_hits <- read_excel("final_filtered_icam_binding.xlsx")

# Vector of samples that had positive ICAM1 hits
icam1_samples <- unique(icam_hits$Sample)

# DC4 samples: motif_alt_id == "ICAM1_1"
dc4_samples <- unique(icam_hits$Sample[icam_hits$motif_alt_id == "ICAM1_1"])

# -------------------------
# Subset expression matrix
# -------------------------
# v$E: voom-normalized expression matrix (rows = genes/domains, cols = samples)
dblb_expr <- v$E[rownames(v$E) %in% icam1_domains, ]

# Keep only samples present in expression matrix AND have CM info AND ICAM1 hits
samples_to_keep <- intersect(colnames(dblb_expr), names(disease_classifications$cerebral))
samples_to_keep <- intersect(samples_to_keep, icam1_samples)
dblb_expr <- dblb_expr[, samples_to_keep]

# -------------------------
# Prepare long-format dataframe
# -------------------------
dblb_df <- as.data.frame(dblb_expr)
dblb_df$GeneID <- rownames(dblb_expr)

dblb_long <- dblb_df %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "expression") %>%
  # Keep only valid samples
  filter(sample %in% samples_to_keep) %>%
  mutate(
    CM_status = ifelse(disease_classifications$cerebral[sample],
                       "severe_cerebral", "severe_nonCerebral"),
    DC4_status = ifelse(sample %in% dc4_samples, "DC4", "Non-DC4")
  ) %>%
  # Drop any NAs just in case
  filter(!is.na(CM_status), !is.na(DC4_status)) %>%
  # Turn into factors with explicit levels
  mutate(
    CM_status = factor(CM_status, levels = c("severe_nonCerebral", "severe_cerebral")),
    DC4_status = factor(DC4_status, levels = c("Non-DC4", "DC4"))
  )

# -------------------------
# Aggregate across domains per sample
# -------------------------
dblb_agg <- dblb_long %>%
  group_by(sample, CM_status, DC4_status) %>%
  summarise(expression_sum = sum(expression), .groups = "drop")

# -------------------------
# Plot boxplot
# -------------------------
plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

ggplot(dblb_agg, aes(x = CM_status, y = expression_sum, fill = DC4_status)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("DC4" = "#d7191c", "Non-DC4" = "#2c7bb6")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) +
  labs(title = "ICAM1-binding DBLb domain expression",
       x = "Cerebral Malaria status", y = "Aggregated logCPM expression")

ggsave(file.path(plots_dir, "ICAM1_DBLb_CM_DC4_boxplot.pdf"), width = 10, height = 6)

####plotting interdomain categories by malaria status####
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# -------------------------
# Parameters and Data Loading
# -------------------------

# Load interdomain mapping
interdomain_mapping <- read_excel("interdomain_mapping.xlsx")

# Create named vector for easy mapping
domain_to_interdomain <- setNames(interdomain_mapping$interdomain_category, 
                                  interdomain_mapping$domain)

# Get unique interdomain categories
interdomain_categories <- unique(interdomain_mapping$interdomain_category)

# -------------------------
# Subset expression matrix
# -------------------------

# Debug: Check what we have
cat("Total domains in v$E:", nrow(v$E), "\n")
cat("Total domains in interdomain mapping:", nrow(interdomain_mapping), "\n")

# Get domains that are in our mapping
domains_to_analyze <- intersect(rownames(v$E), interdomain_mapping$domain)
cat("Domains found with exact matches:", length(domains_to_analyze), "\n")

# Also find domains that match the mapping pattern (e.g., NTSA1 matches NTSA)
additional_domains <- c()
for(mapped_domain in interdomain_mapping$domain) {
  # Find domains in v$E that start with the mapped domain name
  pattern_matches <- rownames(v$E)[grepl(paste0("^", mapped_domain), rownames(v$E))]
  additional_domains <- c(additional_domains, pattern_matches)
}

# Combine exact matches and pattern matches
domains_to_analyze <- unique(c(domains_to_analyze, additional_domains))
cat("Total domains found (exact + pattern matches):", length(domains_to_analyze), "\n")

# Create an updated mapping for pattern matches
updated_domain_to_interdomain <- domain_to_interdomain  # start with exact matches

for(domain in domains_to_analyze) {
  if(!domain %in% names(updated_domain_to_interdomain)) {
    # Find which mapping category this domain belongs to
    for(mapped_domain in interdomain_mapping$domain) {
      if(grepl(paste0("^", mapped_domain), domain)) {
        updated_domain_to_interdomain[domain] <- interdomain_mapping$interdomain_category[interdomain_mapping$domain == mapped_domain]
        break
      }
    }
  }
}

cat("Updated mapping created for", length(updated_domain_to_interdomain), "domains\n")

if(length(domains_to_analyze) == 0) {
  cat("No matching domains found! Checking first few domain names:\n")
  cat("First 10 domains in v$E:", head(rownames(v$E), 10), "\n")
  cat("First 10 domains in mapping:", head(interdomain_mapping$domain, 10), "\n")
  stop("No matching domains found between expression matrix and interdomain mapping")
}

# Subset expression matrix to only these domains
interdomain_expr <- v$E[domains_to_analyze, , drop = FALSE]

# Get samples with disease classification info
# Use the severity factor which has the S/U information
disease_names <- names(disease_classifications$severity)

samples_to_keep <- intersect(colnames(interdomain_expr), disease_names)
cat("Samples with both expression data and disease info:", length(samples_to_keep), "\n")

if(length(samples_to_keep) == 0) {
  cat("No matching samples found!\n")
  cat("First 10 samples in expression matrix:", head(colnames(interdomain_expr), 10), "\n")
  cat("First 10 samples in disease classifications:", head(disease_names, 10), "\n")
  stop("No matching samples found between expression matrix and disease classifications")
}

interdomain_expr <- interdomain_expr[, samples_to_keep, drop = FALSE]

# -------------------------
# Prepare long-format dataframe
# -------------------------

interdomain_df <- as.data.frame(interdomain_expr)
interdomain_df$GeneID <- rownames(interdomain_expr)

# Add interdomain category mapping using updated mapping
interdomain_df$interdomain_category <- updated_domain_to_interdomain[interdomain_df$GeneID]

# Check if we have samples to pivot
cat("Dataframe dimensions before pivot:", nrow(interdomain_df), "x", ncol(interdomain_df), "\n")
cat("Columns to pivot:", ncol(interdomain_df) - 2, "\n")  # minus GeneID and interdomain_category

if(ncol(interdomain_df) <= 2) {
  stop("No sample columns found to pivot. Check sample filtering.")
}

interdomain_long <- interdomain_df %>%
  pivot_longer(-c(GeneID, interdomain_category), names_to = "sample", values_to = "expression") %>%
  # Keep only valid samples
  filter(sample %in% samples_to_keep) %>%
  mutate(
    # Determine malaria status using the severity factor
    malaria_status = case_when(
      disease_classifications$severity[sample] == "Severe" ~ "S",
      disease_classifications$severity[sample] == "Asymptomatic" ~ "U",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Drop any NAs
  filter(!is.na(malaria_status), !is.na(interdomain_category)) %>%
  # Turn into factors
  mutate(
    malaria_status = factor(malaria_status, levels = c("U", "S")),
    interdomain_category = factor(interdomain_category)
  )

# -------------------------
# Debug interdomain categories
# -------------------------

cat("\n=== Debugging interdomain categories ===\n")
cat("All unique interdomain categories in data:\n")
print(unique(interdomain_long$interdomain_category))

cat("\nCategories containing 'NTS':\n")
nts_categories <- unique(interdomain_long$interdomain_category)[grepl("NTS", unique(interdomain_long$interdomain_category))]
print(nts_categories)

cat("\nOriginal domains mapped to NTS categories:\n")
nts_domains <- interdomain_mapping[grepl("NTS", interdomain_mapping$interdomain_category), ]
print(nts_domains)

cat("\nDomains in expression matrix containing 'NTS':\n")
nts_expr_domains <- rownames(v$E)[grepl("NTS", rownames(v$E))]
print(head(nts_expr_domains, 20))  # show first 20

# -------------------------
# Aggregate across domains per sample and interdomain category
# -------------------------

interdomain_agg <- interdomain_long %>%
  group_by(sample, malaria_status, interdomain_category) %>%
  summarise(expression_sum = sum(expression), .groups = "drop")

# -------------------------
# Statistical Tests - Wilcoxon tests within each malaria status and domain family
# -------------------------

# Function to extract domain family from interdomain category
get_domain_family <- function(category) {
  case_when(
    grepl("^CIDR", category) ~ "CIDR",
    grepl("^DBL", category) ~ "DBL",
    grepl("NTS", category) ~ "NTS",  # Changed from ^NTS to just NTS to catch NTSpam
    TRUE ~ "Other"
  )
}

# Function to perform pairwise Wilcoxon tests within domain families
perform_pairwise_wilcox <- function(data, status) {
  cat("\n=== Wilcoxon tests for", status, "malaria (within domain families) ===\n")
  
  # Get data for this status
  status_data <- data %>% 
    filter(malaria_status == status) %>%
    mutate(domain_family = get_domain_family(as.character(interdomain_category)))
  
  # Create results dataframe
  results <- data.frame()
  
  # Get unique domain families
  families <- unique(status_data$domain_family)
  families <- families[families != "Other"]  # exclude "Other" if any
  
  # Perform pairwise comparisons within each domain family
  for(family in families) {
    cat("\n--- ", family, " family comparisons ---\n")
    
    # Get categories within this family
    family_data <- status_data %>% filter(domain_family == family)
    categories <- unique(family_data$interdomain_category)
    
    if(length(categories) < 2) {
      cat("Only", length(categories), "category in", family, "family - no comparisons possible\n")
      next
    }
    
    # Perform pairwise comparisons within this family
    for(i in 1:(length(categories)-1)) {
      for(j in (i+1):length(categories)) {
        cat1 <- categories[i]
        cat2 <- categories[j]
        
        # Get expression values for each category
        expr1 <- family_data %>% filter(interdomain_category == cat1) %>% pull(expression_sum)
        expr2 <- family_data %>% filter(interdomain_category == cat2) %>% pull(expression_sum)
        
        # Perform Wilcoxon test
        if(length(expr1) > 0 & length(expr2) > 0) {
          test_result <- wilcox.test(expr1, expr2)
          
          # Store results
          result_row <- data.frame(
            malaria_status = status,
            domain_family = family,
            category1 = cat1,
            category2 = cat2,
            p_value = test_result$p.value,
            statistic = test_result$statistic
          )
          results <- rbind(results, result_row)
          
          # Print results
          cat(sprintf("%s vs %s: p-value = %.4f, W = %.1f\n", 
                      cat1, cat2, test_result$p.value, test_result$statistic))
        }
      }
    }
  }
  
  return(results)
}

# Perform tests for both severe and asymptomatic
wilcox_results_S <- perform_pairwise_wilcox(interdomain_agg, "S")
wilcox_results_U <- perform_pairwise_wilcox(interdomain_agg, "U")

# Combine results
all_wilcox_results <- rbind(wilcox_results_S, wilcox_results_U)

# Apply multiple testing correction
all_wilcox_results$p_adjusted <- p.adjust(all_wilcox_results$p_value, method = "BH")

# Print adjusted results
cat("\n=== Results with multiple testing correction (BH) ===\n")
cat("Results grouped by domain family:\n")
for(family in unique(all_wilcox_results$domain_family)) {
  cat("\n", family, "family:\n")
  family_results <- all_wilcox_results %>% filter(domain_family == family)
  print(family_results[, c("malaria_status", "category1", "category2", "p_value", "p_adjusted")])
}

# -------------------------
# Create plots directory
# -------------------------
plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

# -------------------------
# Boxplot
# -------------------------

p1 <- ggplot(interdomain_agg, aes(x = interdomain_category, y = expression_sum, fill = malaria_status)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("S" = "#d7191c", "U" = "#2c7bb6"),
                    labels = c("S" = "Severe", "U" = "Asymptomatic")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Interdomain category expression by malaria status",
       x = "Interdomain category", 
       y = "Aggregated logCPM expression",
       fill = "Malaria Status")

ggsave(file.path(plots_dir, "interdomain_malaria_status_boxplot.pdf"), 
       plot = p1, width = 12, height = 8)

# -------------------------
# Separate plots for each malaria status
# -------------------------

# Severe malaria only
p2 <- interdomain_agg %>%
  filter(malaria_status == "S") %>%
  ggplot(aes(x = interdomain_category, y = expression_sum)) +
  geom_boxplot(fill = "#d7191c", alpha = 0.7, outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Interdomain category expression - Severe Malaria",
       x = "Interdomain category", 
       y = "Aggregated logCPM expression")

ggsave(file.path(plots_dir, "interdomain_severe_malaria_boxplot.pdf"), 
       plot = p2, width = 10, height = 6)

# Asymptomatic malaria only
p3 <- interdomain_agg %>%
  filter(malaria_status == "U") %>%
  ggplot(aes(x = interdomain_category, y = expression_sum)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.7, outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Interdomain category expression - Asymptomatic Malaria",
       x = "Interdomain category", 
       y = "Aggregated logCPM expression")

ggsave(file.path(plots_dir, "interdomain_Asymptomatic_malaria_boxplot.pdf"), 
       plot = p3, width = 10, height = 6)

# -------------------------
# Save statistical results
# -------------------------

# Save Wilcoxon test results
write.csv(all_wilcox_results, file.path(plots_dir, "wilcoxon_test_results.csv"), 
          row.names = FALSE)

# Print summary
cat("\n=== Summary ===\n")
cat("Total samples analyzed:", length(unique(interdomain_agg$sample)), "\n")
cat("Severe malaria samples:", length(unique(interdomain_agg$sample[interdomain_agg$malaria_status == "S"])), "\n")
cat("Asymptomatic malaria samples:", length(unique(interdomain_agg$sample[interdomain_agg$malaria_status == "U"])), "\n")
cat("Interdomain categories:", paste(levels(interdomain_agg$interdomain_category), collapse = ", "), "\n")
cat("Significant comparisons (p < 0.05):", sum(all_wilcox_results$p_adjusted < 0.05), "\n")
###################################################################################################################
# DOMAIN ANALYSIS
###################################################################################################################

domains_SM <- c("CIDRa1_1","CIDRa1_4","CIDRa1_5","CIDRa1_8","CIDRa3_1","CIDRa3_4","CIDRb2","CIDRb5","DBLa0_1","DBLa0_16",
                "DBLa1_2","DBLa1_7","DBLa2","DBLb12","DBLb3","DBLb5","DBLb6","DBLb7","DBLd1","DBLd5",
                "DBLe11","DBLe4","DBLe6","DBLg11","DBLg12","DBLg2","DBLg4","DBLg5","DBLg6","DBLg7","DBLpam1","DBLz4","DBLz5","NTSB3","NTSA5")

domains_UM <- c("CIDRa1_2","CIDRa1_3","CIDRa2_11","CIDRa2_5","CIDRa3_3",
                "CIDRd2","CIDRg12","CIDRg8","DBLa0_15","DBLa0_23",
                "DBLa0_24","DBLb13","DBLb2","DBLb9","DBLd8","DBLd9","DBLe12","DBLe14",
                "DBLg16","DBLg18","DBLz1","NTSA1", "NTSB6")

cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
severe_samples <- colnames(cpm_matrix)[grepl("^CH|^ch", colnames(cpm_matrix))]
Asymptomatic_samples <- setdiff(colnames(cpm_matrix), severe_samples)

min_cpm <- 1; min_samples <- 3; min_max <- 5
check_expr <- function(dom, samples){
  if(!dom %in% rownames(cpm_matrix)) return(FALSE)
  vals <- cpm_matrix[dom, samples]
  sum(vals >= min_cpm) >= min_samples & max(vals) >= min_max
}

filtered_SM <- domains_SM[sapply(domains_SM, check_expr, samples = severe_samples)]
filtered_UM <- domains_UM[sapply(domains_UM, check_expr, samples = Asymptomatic_samples)]

find_top_hits <- function(domains, samples, top_n = 10, min_cpm_req = 1){
  top_hits <- lapply(domains, function(dom){
    if(!dom %in% rownames(cpm_matrix)) return(NULL)
    vals <- cpm_matrix[dom, samples]
    vals <- vals[vals >= min_cpm_req]
    if(length(vals)==0) return(NULL)
    n_take <- min(top_n, length(vals))
    data.frame(Domain=dom, Sample=names(sort(vals, decreasing=TRUE)[1:n_take]),
               CPM=sort(vals, decreasing=TRUE)[1:n_take], stringsAsFactors=FALSE)
  })
  do.call(rbind, top_hits[sapply(top_hits, Negate(is.null))])
}

top_SM <- find_top_hits(filtered_SM, severe_samples)
top_UM <- find_top_hits(filtered_UM, Asymptomatic_samples)

create_mapping <- function(top_hits){
  sapply(unique(top_hits$Sample), function(s) top_hits$Domain[top_hits$Sample == s], simplify=FALSE)
}

greedy_set <- function(mapping, domains){
  covered <- c(); selected <- c(); pool <- mapping
  while(length(covered) < length(intersect(domains, unique(unlist(mapping)))) && length(pool)>0){
    scores <- sapply(pool, function(x) length(intersect(setdiff(x, covered), domains)))
    if(max(scores)==0) break
    best <- names(which.max(scores))
    selected <- c(selected, best)
    covered <- union(covered, pool[[best]])
    pool[[best]] <- NULL
  }
  selected
}

map_SM <- create_mapping(top_SM); map_UM <- create_mapping(top_UM)
min_SM <- greedy_set(map_SM, filtered_SM); min_UM <- greedy_set(map_UM, filtered_UM)

assign_domains <- function(domains, min_samples, top_hits){
  do.call(rbind, lapply(domains, function(dom){
    candidates <- top_hits[top_hits$Domain==dom & top_hits$Sample %in% min_samples, ]
    if(nrow(candidates)==0) return(NULL)
    best <- candidates[which.max(candidates$CPM), ]
    data.frame(Domain=dom, Primary_Sample=best$Sample, CPM=best$CPM, stringsAsFactors=FALSE)
  }))
}

assign_SM <- assign_domains(filtered_SM, min_SM, top_SM)
assign_UM <- assign_domains(filtered_UM, min_UM, top_UM)
assign_all <- rbind(assign_SM, assign_UM)
write.csv(assign_all, file.path(results_dir, "domain_sample_assignments.csv"), row.names=FALSE)

summary_stats <- data.frame(
  Group=c("Severe","Asymptomatic","Combined"),
  Domains=c(length(filtered_SM), length(filtered_UM), length(filtered_SM)+length(filtered_UM)),
  Samples_Required=c(length(min_SM), length(min_UM), length(min_SM)+length(min_UM)),
  Total_Samples=c(length(severe_samples), length(Asymptomatic_samples), length(severe_samples)+length(Asymptomatic_samples))
)
write.csv(summary_stats, file.path(results_dir, "severity_analysis_summary.csv"), row.names=FALSE)

###################################################################################################################
# RLE PLOTS (Relative Log Expression)
###################################################################################################################

plotVoomRLE <- function(E, colours, title = "RLE Plot"){
  mn <- apply(E, 1, median) 
  rle <- data.frame(sweep(E, MARGIN=1, STATS=mn, FUN='-')) 
  boxplot(rle, col=colours, outline=FALSE, las=2, ylim=c(-7,7), main=title)
  abline(h = 0, col = "black")
}

# RLE plot for normalized data
pdf(file.path(plots_dir, "RLE_plot_normalized.pdf"), width = 16, height = 8)
plotVoomRLE(v$E, colors[categories], "RLE Plot - Normalized Data")
dev.off()

###################################################################################################################
# Additional PCA plots to determine clustering within syndromes - UPDATED WITH REAL RUV4
###################################################################################################################

disease_codes <- dge$samples$Disease
new_classifications <- rep("Other", length(disease_codes))
names(new_classifications) <- rownames(dge$samples)

new_classifications[grepl("C", disease_codes)] <- "Cerebral"
acidosis_hypo <- (grepl("L", disease_codes) | grepl("H", disease_codes)) & !grepl("C", disease_codes)
new_classifications[acidosis_hypo] <- "Acidosis_Hypo"
anaemia <- grepl("A", disease_codes) & !grepl("C", disease_codes)
new_classifications[anaemia] <- "Anaemia"
new_classifications[disease_codes == "U"] <- "Asymptomatic"

grp_detailed <- factor(new_classifications,
                       levels = c("Asymptomatic", "Cerebral", "Acidosis_Hypo", "Anaemia", "Other"))

colors_detailed <- c("Asymptomatic" = "#d7191c", "Cerebral" = "#008080", 
                     "Acidosis_Hypo" = "#FFA500", "Anaemia" = "#A020F0", 
                     "Other" = "#808080")

# Stage covariates for detailed PCA
stage_covs_list_detailed <- list(
  Ring = covs_base$Ring,
  Gametocyte = covs_base$Gametocyte,
  Schizont = covs_base$Schizont
)

# Use RUV4 with specified number of factors
k_factors <- 2  # Change this value based on supervisor's decision (tested k=1,2,3 - see k_factor_comparison_summary.csv)

cat("Creating PCA plots for k =", k_factors, "RUV factors\n")

# Use severity grouping instead of detailed grouping to avoid singularity
grp_simple <- disease_classifications$severity[colnames(v$E)]
grp_simple <- grp_simple[!is.na(grp_simple)]
expr_filtered <- v$E[, names(grp_simple)]

# Use real RUV4 with specified k value on filtered data
ruv_factors_detailed <- compute_real_ruv4(expr_filtered, grp_simple, empirical_controls_genes, k_factors = k_factors)
corrected_expr_detailed <- get_corrected_expression(expr_filtered, ruv_factors_detailed$W, stage_covs_list_detailed)

# PCA on all genes
pca_detailed <- autoplot(
  prcomp(t(corrected_expr_detailed)),
  data = data.frame(Group = grp_detailed[colnames(corrected_expr_detailed)]),
  colour = 'Group',
  size = 1,
  alpha = 0.5
) +
  scale_color_manual(values = colors_detailed) +
  ggtitle(paste("PCA: Library Size + All Stages + RUV (k =", k_factors, ") - Detailed Disease Classification")) +
  theme_bw()

ggsave(file.path(plots_dir, "PCA_detailed_disease_classification.pdf"), pca_detailed, width = 10, height = 8)

# Identify var genes (not starting with PF3D7)
var_genes <- !grepl("^PF3D7", rownames(expr_filtered))

# Subset expression matrix to only var genes
expr_var <- expr_filtered[var_genes, ]

# Recompute RUV factors on var-only data
ruv_factors_var <- compute_real_ruv4(expr_var, grp_simple, empirical_controls_genes, k_factors = k_factors)
corrected_expr_var <- get_corrected_expression(expr_var, ruv_factors_var$W, stage_covs_list_detailed)

# PCA on var-only data
pca_var <- autoplot(
  prcomp(t(corrected_expr_var)),
  data = data.frame(Group = grp_detailed[colnames(corrected_expr_var)]),
  colour = 'Group',
  size = 1,
  alpha = 0.5
) +
  scale_color_manual(values = colors_detailed) +
  ggtitle("PCA: Var Gene Expression (Detailed Disease Classification)") +
  theme_bw()

ggsave(file.path(plots_dir, "PCA_var_only_disease_classification.pdf"), pca_var, width = 10, height = 8)

# Optionally export var-only matrices
write.csv(corrected_expr_var, file.path(results_dir, "corrected_var_expression.csv"))
write.csv(expr_var, file.path(results_dir, "uncorrected_var_log2cpm.csv"))

# Export corrected and uncorrected expression data for individual samples
write.csv(corrected_expr_detailed, 
          file.path(results_dir, "corrected_expression_individual_samples.csv"))
write.csv(v$E, file.path(results_dir, "uncorrected_log2cpm_expression.csv"))

# Also export sample metadata
write.csv(dge$samples, 
          file.path(results_dir, "sample_metadata.csv"))

###################################################################################################################
# Comparing number of factors in RUV4 analysis
###################################################################################################################

library(ruv)

compute_real_ruv4 <- function(expression_data, group_factor, control_genes, k_factors = 3) {
  
  # Convert group factor to numeric matrix (like original pipeline)
  categoriesRUV <- data.matrix(as.numeric(group_factor == "Disease"))
  if(length(unique(group_factor)) == 2 && "Severe" %in% levels(group_factor)) {
    categoriesRUV <- data.matrix(as.numeric(group_factor == "Severe"))
  }
  
  # Create empirical controls logical vector
  empirical_controls <- rownames(expression_data) %in% control_genes
  
  # Transpose expression data (genes become rows x samples columns -> samples x genes)
  genes <- data.matrix(t(expression_data))
  
  # Run RUV4 (no Z matrix for now - add staging separately)
  ruv_result <- RUV4(genes, categoriesRUV, empirical_controls, k_factors)
  
  return(list(W = ruv_result$W))
}

# MODIFIED ANALYSIS LOOP WITH K-FACTOR COMPARISON
all_results <- list(); all_plots <- list()

# Test different k values
k_values <- c(1, 2, 3)

for(disease_name in names(disease_classifications)) {
  selected_class <- disease_classifications[[disease_name]]
  selected_class_vector <- as.vector(unlist(selected_class))
  
  if(is.logical(selected_class_vector)) {
    valid <- !is.na(selected_class_vector)
    filtered_expression <- v$E[, valid, drop=FALSE]
    filtered_covs <- covs_base[valid, , drop=FALSE]
    grp <- factor(ifelse(selected_class_vector[valid], "Disease", "Asymptomatic"),
                  levels = c("Asymptomatic", "Disease"))
  } else if(disease_name == "severity") {
    valid <- !is.na(selected_class_vector)
    filtered_expression <- v$E[, valid, drop=FALSE]
    filtered_covs <- covs_base[valid, , drop=FALSE]
    grp <- selected_class[valid]
  }
  
  # Stage covariates (same as before)
  covs_use <- filtered_covs
  for(stage_col in c("Ring", "Gametocyte", "Schizont")) {
    if(stage_col %in% colnames(covs_use) && is.numeric(covs_use[[stage_col]])) {
      covs_use[[stage_col]] <- scale(covs_use[[stage_col]], center = TRUE, scale = FALSE)
    }
  }
  
  df_base <- data.frame(disease = grp, covs_use)
  D0 <- model.matrix(~ 0 + disease, data = df_base)
  D_ring <- cbind(D0, Ring = df_base$Ring)
  D_ring_all <- cbind(D0, Ring = df_base$Ring, Gametocyte = df_base$Gametocyte, Schizont = df_base$Schizont)
  
  colnames(D0) <- make.names(colnames(D0))
  colnames(D_ring) <- make.names(colnames(D_ring))
  colnames(D_ring_all) <- make.names(colnames(D_ring_all))
  
  # Test different k values for RUV
  for(k in k_values) {
    cat("Computing RUV factors for", disease_name, "with k =", k, "...\n")
    ruv_factors <- compute_real_ruv4(filtered_expression, grp, empirical_controls_genes, k_factors = k)
    
    D_ruv_no <- cbind(D0, ruv_factors$W)
    D_ring_ruv <- cbind(D_ring, ruv_factors$W)
    D_all_ruv <- cbind(D_ring_all, ruv_factors$W)
    
    design_list <- list(D0, D_ruv_no, D_ring, D_ring_ruv, D_ring_all, D_all_ruv)
    method_names <- c("Library_Size_Only", paste0("Library_Size_RUV_k", k), "Library_Size_Ring", 
                      paste0("Library_Size_Ring_RUV_k", k), "Library_Size_All_Stages", 
                      paste0("Library_Size_All_Stages_RUV_k", k))
    
    limma_results_list <- list(); volcano_plots <- list(); pca_plots <- list(); ma_plots <- list()
    
    for(i in seq_along(design_list)) {
      design <- design_list[[i]]
      fit <- limma::lmFit(filtered_expression, design)
      
      design_cols <- colnames(design)
      if(disease_name == "severity") {
        disease_col_fixed <- design_cols[grepl("Severe", design_cols)][1]
        control_col_fixed <- design_cols[grepl("Asymptomatic", design_cols)][1]
      } else {
        disease_col_fixed <- design_cols[grepl("Disease", design_cols)][1]
        control_col_fixed <- design_cols[grepl("Asymptomatic", design_cols)][1]
      }
      
      if(!is.na(disease_col_fixed) && !is.na(control_col_fixed)) {
        contrast_matrix <- matrix(0, nrow = ncol(design), ncol = 1)
        rownames(contrast_matrix) <- colnames(design)
        colnames(contrast_matrix) <- "Disease_vs_Control"
        contrast_matrix[disease_col_fixed, 1] <- 1
        contrast_matrix[control_col_fixed, 1] <- -1
        fit2 <- limma::contrasts.fit(fit, contrast_matrix)
        fit2 <- limma::eBayes(fit2)
        results <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
      } else {
        fit2 <- limma::eBayes(fit)
        coef_idx <- which(grepl("Disease|Severe", colnames(fit2$coefficients)))[1]
        results <- limma::topTable(fit2, coef = ifelse(!is.na(coef_idx), coef_idx, 2), 
                                   number = Inf, sort.by = "none")
      }
      
      results$gene <- rownames(results)
      limma_results_list[[i]] <- results
      
      # Create plots
      volcano_title <- paste("Volcano:", gsub("_", " ", method_names[i]), "-", disease_name)
      volcano_plots[[i]] <- create_volcano_plot(results, volcano_title, FC_THRESHOLD, PVAL_THRESHOLD)
      
      # Corrected expression for visualization
      stage_covs_list <- NULL; ruv_use <- NULL
      if (grepl("Ring", method_names[i]) && !grepl("All_Stages", method_names[i])) {
        stage_covs_list <- list(Ring = df_base$Ring)
      } else if (grepl("All_Stages", method_names[i])) {
        stage_covs_list <- list(Ring = df_base$Ring, Gametocyte = df_base$Gametocyte, Schizont = df_base$Schizont)
      }
      if (grepl("RUV", method_names[i])) ruv_use <- ruv_factors$W
      
      corrected_expr <- get_corrected_expression(filtered_expression, ruv_use, stage_covs_list)
      pca_plots[[i]] <- create_pca_plot(corrected_expr, grp, paste("PCA:", gsub("_", " ", method_names[i])))
      ma_plots[[i]] <- create_ma_plot(results, paste("MA:", gsub("_", " ", method_names[i])))
      
      # Save results
      write.csv(results, file.path(results_dir, paste0(method_names[i], "_", disease_name, "_results.csv")))
    }
    
    # Save plots for this k value
    pdf(file.path(plots_dir, paste0("PCA_grid_", disease_name, "_k", k, ".pdf")), width = 12, height = 12)
    grid.arrange(grobs = pca_plots, ncol = 2, nrow = 3, top = paste("PCA Comparison:", disease_name, "k =", k))
    dev.off()
    
    pdf(file.path(plots_dir, paste0("Volcano_grid_", disease_name, "_k", k, ".pdf")), width = 12, height = 12)
    grid.arrange(grobs = volcano_plots, ncol = 2, nrow = 3, top = paste("Volcano Comparison:", disease_name, "k =", k))
    dev.off()
    
    # Store results for this k value
    all_results[[paste(disease_name, "k", k, sep = "_")]] <- limma_results_list
    all_plots[[paste(disease_name, "k", k, sep = "_")]] <- volcano_plots
  }
}

# CREATE K-FACTOR COMPARISON SUMMARY
k_comparison_summary <- data.frame()

for(disease_name in names(disease_classifications)) {
  for(k in k_values) {
    key <- paste(disease_name, "k", k, sep = "_")
    if(key %in% names(all_results)) {
      # Get results from Library_Size_All_Stages_RUV (method 6)
      limma_results <- all_results[[key]][[6]]
      if(!is.null(limma_results)) {
        sig_genes <- sum(!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1)
        core_sig <- sum(!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1 & 
                          grepl("^PF3D7", limma_results$gene))
        var_sig <- sum(!is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.1 & 
                         !grepl("^PF3D7", limma_results$gene))
        
        summary_row <- data.frame(
          Disease = disease_name,
          K_Factors = k,
          Total_Significant = sig_genes,
          Core_Significant = core_sig,
          Var_Significant = var_sig
        )
        k_comparison_summary <- rbind(k_comparison_summary, summary_row)
      }
    }
  }
}

write.csv(k_comparison_summary, file.path(results_dir, "k_factor_comparison_summary.csv"), row.names = FALSE)
print(k_comparison_summary)

cat("K-factor comparison completed!\n")
cat("Check the summary table to see how different k values affect significant gene detection.\n")

##########################################
cat("Analysis completed successfully!\n")
cat("Results saved to:", results_dir, "\n")
cat("Plots saved to:", plots_dir, "\n")
cat("GSEA files saved to:", gsea_dir, "\n")