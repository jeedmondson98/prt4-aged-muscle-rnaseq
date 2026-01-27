#!/usr/bin/env Rscript

# =============================================================================
# tximport Script for Salmon Quantification Data
# Aggregates transcript-level counts to gene-level with gene symbols
# =============================================================================

# Load required packages (install if needed)
required_packages <- c("tximport", "readr", "dplyr", "biomaRt", "GenomicFeatures")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("tximport", "biomaRt", "GenomicFeatures")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# Configuration
# =============================================================================

# Base directory for Salmon output
salmon_base <- "/Volumes/je_toshiba/salmon_output"

# Output directory for gene-level counts
output_dir <- "/Volumes/je_toshiba/gene_counts"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Create transcript-to-gene mapping using biomaRt
# =============================================================================

message("Creating transcript-to-gene mapping from Ensembl...")

# Connect to Ensembl for rat
ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Get transcript to gene mapping with gene symbols
tx2gene <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
  mart = ensembl
)

# Rename columns
colnames(tx2gene) <- c("TXNAME", "GENEID", "SYMBOL")

# Handle missing gene symbols - use gene ID if no symbol
tx2gene$SYMBOL[tx2gene$SYMBOL == ""] <- tx2gene$GENEID[tx2gene$SYMBOL == ""]

# Save tx2gene mapping for reference
write.csv(tx2gene, file.path(output_dir, "tx2gene_mapping.csv"), row.names = FALSE)
message(paste("Saved transcript-to-gene mapping:", nrow(tx2gene), "transcripts"))

# =============================================================================
# Function to process a dataset
# =============================================================================

process_dataset <- function(dataset_name, salmon_dir, tx2gene, output_dir) {
  
  message(paste("\n=== Processing", dataset_name, "==="))
  
  # Find all sample directories
  sample_dirs <- list.dirs(salmon_dir, full.names = TRUE, recursive = FALSE)
  sample_names <- basename(sample_dirs)
  
  message(paste("Found", length(sample_names), "samples"))
  
  # Build file paths to quant.sf files
  files <- file.path(sample_dirs, "quant.sf")
  names(files) <- sample_names
  
  # Check all files exist
  missing <- !file.exists(files)
  if (any(missing)) {
    warning(paste("Missing quant.sf files:", paste(names(files)[missing], collapse = ", ")))
    files <- files[!missing]
  }
  
  # Use countsFromAbundance = "no" for limma-voom
  # voom handles normalization internally, so raw estimated counts are preferred
  txi <- tximport(files, 
                  type = "salmon", 
                  tx2gene = tx2gene[, c("TXNAME", "GENEID")],
                  ignoreTxVersion = TRUE,
                  countsFromAbundance = "no")
  
  message(paste("Imported", nrow(txi$counts), "genes across", ncol(txi$counts), "samples"))
  
  # =============================================================================
  # Create output matrices with gene symbols
  # =============================================================================
  
  # Get unique gene ID to symbol mapping (use dplyr:: to avoid AnnotationDbi conflict)
  gene_symbols <- tx2gene %>%
    dplyr::select(GENEID, SYMBOL) %>%
    dplyr::distinct(GENEID, .keep_all = TRUE)
  
  # Create counts matrix with gene symbols
  counts_df <- as.data.frame(txi$counts)
  counts_df$ensembl_gene_id <- rownames(counts_df)
  counts_df <- merge(counts_df, gene_symbols, by.x = "ensembl_gene_id", by.y = "GENEID", all.x = TRUE)
  
  # Handle any missing symbols
  counts_df$SYMBOL[is.na(counts_df$SYMBOL)] <- counts_df$ensembl_gene_id[is.na(counts_df$SYMBOL)]
  
  # Reorder columns: gene_symbol, ensembl_id, then samples
  sample_cols <- setdiff(colnames(counts_df), c("ensembl_gene_id", "SYMBOL"))
  counts_df <- counts_df[, c("SYMBOL", "ensembl_gene_id", sample_cols)]
  colnames(counts_df)[1:2] <- c("gene_symbol", "ensembl_gene_id")
  
  # Create TPM matrix with gene symbols
  tpm_df <- as.data.frame(txi$abundance)
  tpm_df$ensembl_gene_id <- rownames(tpm_df)
  tpm_df <- merge(tpm_df, gene_symbols, by.x = "ensembl_gene_id", by.y = "GENEID", all.x = TRUE)
  tpm_df$SYMBOL[is.na(tpm_df$SYMBOL)] <- tpm_df$ensembl_gene_id[is.na(tpm_df$SYMBOL)]
  tpm_df <- tpm_df[, c("SYMBOL", "ensembl_gene_id", sample_cols)]
  colnames(tpm_df)[1:2] <- c("gene_symbol", "ensembl_gene_id")
  
  # Create length matrix with gene symbols
  length_df <- as.data.frame(txi$length)
  length_df$ensembl_gene_id <- rownames(length_df)
  length_df <- merge(length_df, gene_symbols, by.x = "ensembl_gene_id", by.y = "GENEID", all.x = TRUE)
  length_df$SYMBOL[is.na(length_df$SYMBOL)] <- length_df$ensembl_gene_id[is.na(length_df$SYMBOL)]
  length_df <- length_df[, c("SYMBOL", "ensembl_gene_id", sample_cols)]
  colnames(length_df)[1:2] <- c("gene_symbol", "ensembl_gene_id")
  
  # =============================================================================
  # Save outputs
  # =============================================================================
  
  # Create dataset output directory
  dataset_out <- file.path(output_dir, dataset_name)
  dir.create(dataset_out, showWarnings = FALSE)
  
  # Save as TSV (tab-separated, easy to read)
  write.table(counts_df, file.path(dataset_out, "gene_counts.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(tpm_df, file.path(dataset_out, "gene_tpm.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(length_df, file.path(dataset_out, "gene_lengths.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Also save as CSV
  write.csv(counts_df, file.path(dataset_out, "gene_counts.csv"), row.names = FALSE)
  write.csv(tpm_df, file.path(dataset_out, "gene_tpm.csv"), row.names = FALSE)
  
  # Save the tximport object for downstream analysis (DESeq2, limma-voom)
  saveRDS(txi, file.path(dataset_out, "txi_object.rds"))
  
  message(paste("Saved outputs to:", dataset_out))
  message(paste("  - gene_counts.txt:", nrow(counts_df), "genes x", length(sample_cols), "samples"))
  message(paste("  - gene_tpm.txt (normalized)"))
  message(paste("  - gene_lengths.txt (for offset correction)"))
  message(paste("  - txi_object.rds (for DESeq2/limma)"))
  
  return(txi)
}

# =============================================================================
# Process both datasets
# =============================================================================

# Process PRT2_4 (your main dataset)
txi_prt2_4 <- process_dataset(
  dataset_name = "PRT2_4",
  salmon_dir = file.path(salmon_base, "PRT2_4"),
  tx2gene = tx2gene,
  output_dir = output_dir
)

# Process Shavlakadze2019 (validation dataset)
txi_shavlakadze <- process_dataset(
  dataset_name = "Shavlakadze2019",
  salmon_dir = file.path(salmon_base, "Shavlakadze2019"),
  tx2gene = tx2gene,
  output_dir = output_dir
)

# =============================================================================
# Summary
# =============================================================================

message("\n=============================================================================")
message("TXIMPORT COMPLETE")
message("=============================================================================")
message(paste("Output directory:", output_dir))
message("\nFiles created for each dataset:")
message("  - gene_counts.txt : Raw estimated counts (for limma-voom)")
message("  - gene_tpm.txt    : TPM normalized (for visualization)")
message("  - gene_lengths.txt: Gene lengths (for offset correction if needed)")
message("  - txi_object.rds  : R object (for direct use in DESeq2/limma)")
message("  - gene_counts.csv : CSV format of counts")
message("  - gene_tpm.csv    : CSV format of TPM")
message("\nColumn format: gene_symbol, ensembl_gene_id, sample1, sample2, ...")
message("\nNote: Using countsFromAbundance = 'no' (raw estimated counts)")
message("      This is recommended for limma-voom which handles normalization internally")
message("=============================================================================")


# =============================================================================
# Merge datasets into master file
# =============================================================================

# Read the two count files
prt2_4 <- read.csv("/Volumes/je_toshiba/gene_counts/PRT2_4/gene_counts.csv")
shav <- read.csv("/Volumes/je_toshiba/gene_counts/Shavlakadze2019/gene_counts.csv")

# Merge by gene_symbol and ensembl_gene_id
master_counts <- merge(prt2_4, shav, by = c("gene_symbol", "ensembl_gene_id"), all = TRUE)

# Check dimensions
message(paste("Master file:", nrow(master_counts), "genes x", ncol(master_counts) - 2, "samples"))

# Save
write.csv(master_counts, "/Volumes/je_toshiba/gene_counts/master_gene_counts.csv", row.names = FALSE)
write.table(master_counts, "/Volumes/je_toshiba/gene_counts/master_gene_counts.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Saved master_gene_counts.csv and master_gene_counts.txt")

# =============================================================================
# Create metadata
# =============================================================================

# Get sample names from count files
prt2_4_samples <- setdiff(colnames(prt2_4), c("gene_symbol", "ensembl_gene_id"))
shav_samples <- setdiff(colnames(shav), c("gene_symbol", "ensembl_gene_id"))

# Create animal metadata lookup from your data
animal_meta <- data.frame(
  animal_id = c(421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436,
                437, 440, 441, 442, 443, 444, 455, 453, 445, 446, 448, 449, 450, 451, 452, 454,
                257, 258, 259, 260, 261, 262, 263, 264, 265, 267, 268, 275),
  timepoint = c(20, 20, 20, 20, 60, 60, 60, 60, 30, 30, 30, 30, 30, 30, 30, 30,
                10, 10, 10, 10, 10, 10, 10, 10, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 10, 10, 10, 20, 20, 20, 30, 30, 30),
  sex = c("f", "f", "m", "m", "f", "f", "m", "m", "m", "m", "f", "f", "m", "m", "f", "f",
          "m", "m", "m", "f", "m", "f", "f", "m", "m", "m", "f", "m", "m", "f", "f", "f",
          "m", "m", "m", "m", "m", "m", "m", "m", "m", "m", "m", "m"),
  age_group = c(rep("Old_14mo", 32), rep("Young_3mo", 12)),
  stim = c("Spillover", "Spillover", "Spillover", "Spillover", "CLF", "CLF", "CLF", "CLF",
           "Spillover", "Spillover", "Spillover", "Spillover", "CLF", "CLF", "CLF", "CLF",
           "Spillover", "CLF", "CLF", "Spillover", "CLF", "CLF", "Spillover", "Spillover",
           "Spillover", "Spillover", "Spillover", "CLF", "CLF", "CLF", "CLF", "Spillover",
           rep("Spillover", 12)),
  stringsAsFactors = FALSE
)

# Create PRT2_4 sample metadata
prt2_4_meta <- data.frame(
  sample_id = prt2_4_samples,
  stringsAsFactors = FALSE
)

# Parse sample components
prt2_4_meta$muscle <- ifelse(grepl("^EDL", prt2_4_meta$sample_id), "EDL", "TA")
prt2_4_meta$limb <- ifelse(grepl("L[0-9]", prt2_4_meta$sample_id), "Left", "Right")
prt2_4_meta$animal_id <- as.numeric(gsub(".*([0-9]{3})$", "\\1", prt2_4_meta$sample_id))

# Merge with animal metadata
prt2_4_meta <- merge(prt2_4_meta, animal_meta, by = "animal_id", all.x = TRUE)

# Add dataset info - PRT2 for young, PRT4 for old
prt2_4_meta$dataset <- ifelse(prt2_4_meta$age_group == "Young_3mo", "PRT2", "PRT4")

# EDL samples are from young animals so they're PRT2
prt2_4_meta$dataset[prt2_4_meta$muscle == "EDL"] <- "PRT2"

prt2_4_meta$seq_batch <- ifelse(prt2_4_meta$age_group == "Young_3mo", "2023", "2025")
prt2_4_meta$library_type <- ifelse(prt2_4_meta$age_group == "Young_3mo" & prt2_4_meta$muscle == "TA", "ISR", "IU")

# EDL samples are from young animals but sequenced in 2025
prt2_4_meta$seq_batch[prt2_4_meta$muscle == "EDL"] <- "2025"

# Reorder columns
prt2_4_meta <- prt2_4_meta[, c("sample_id", "animal_id", "muscle", "limb", "age_group", 
                               "timepoint", "sex", "stim", "dataset", "seq_batch", "library_type")]

# Load Shavlakadze metadata from SRA (you may need to adjust this path)
# For now, create placeholder metadata based on sample names
shav_meta_clean <- data.frame(
  sample_id = shav_samples,
  animal_id = seq_along(shav_samples),
  muscle = "Gastrocnemius",
  limb = NA,
  age_group = NA,  # Will need to be filled from SRA metadata
  timepoint = NA,
  sex = NA,
  stim = "None",
  dataset = "Shavlakadze2019",
  seq_batch = "2019",
  library_type = "ISR",
  stringsAsFactors = FALSE
)

# Combine all metadata
master_meta <- rbind(prt2_4_meta, shav_meta_clean)

# Check distributions
message("\n=== Dataset by Age Group ===")
print(table(master_meta$dataset, master_meta$age_group, useNA = "ifany"))

message("\n=== Dataset by Library Type ===")
print(table(master_meta$dataset, master_meta$library_type))

message("\n=== Dataset by Timepoint ===")
print(table(master_meta$dataset, master_meta$timepoint, useNA = "ifany"))

# Save
write.csv(master_meta, "/Volumes/je_toshiba/gene_counts/master_metadata.csv", row.names = FALSE)
write.table(master_meta, "/Volumes/je_toshiba/gene_counts/master_metadata.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\n=== SAVED ===")
message("master_metadata.csv: ", nrow(master_meta), " samples")
message("Location: /Volumes/je_toshiba/gene_counts/")

# =============================================================================
# Copy to project directory
# =============================================================================

# Source directory
source_dir <- "/Volumes/je_toshiba/gene_counts"

# Destination directory
dest_dir <- "/Users/jackedmondson/Library/CloudStorage/OneDrive-LiverpoolJohnMooresUniversity/PhD/Projects/PRT4/data/PRT4_paper_data"

# Create destination if it doesn't exist
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

# Files to copy
files_to_copy <- c(
  "master_gene_counts.csv",
  "master_gene_counts.txt",
  "master_metadata.csv",
  "master_metadata.txt",
  "tx2gene_mapping.csv",
  "PRT2_4/txi_object.rds",
  "PRT2_4/gene_counts.csv",
  "PRT2_4/gene_tpm.csv",
  "Shavlakadze2019/txi_object.rds",
  "Shavlakadze2019/gene_counts.csv",
  "Shavlakadze2019/gene_tpm.csv"
)

# Create subdirectories
dir.create(file.path(dest_dir, "PRT2_4"), showWarnings = FALSE)
dir.create(file.path(dest_dir, "Shavlakadze2019"), showWarnings = FALSE)

# Copy files
for (f in files_to_copy) {
  src <- file.path(source_dir, f)
  dst <- file.path(dest_dir, f)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    message("Copied: ", f)
  } else {
    message("Missing: ", f)
  }
}

# Verify
message("\n=== Files in destination ===")
list.files(dest_dir, recursive = TRUE)