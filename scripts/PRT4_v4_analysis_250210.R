# ================================================================
# PRT4 COMPLETE ANALYSIS SCRIPT v5 - WITH GENE FILTERING
# ================================================================

rm(list = ls())

# ================================================================
# PACKAGES
# ================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readxl)
library(edgeR)
library(limma)
library(pheatmap)
library(mixOmics)
library(fgsea)
library(msigdbr)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(reshape2)
library(dorothea)
library(viper)
library(scales)
library(GSVA)
library(grid)
library(png)
library(cowplot)
library(VennDiagram)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggh4x)
library(ggrepel)

# ================================================================
# SETUP
# ================================================================

setwd("..") # Set to project root

# Create output directories
dir.create("figures/DEGs", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/GSEA", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/sPLSDA", showWarnings = FALSE, recursive = TRUE)
dir.create("paper_v4/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("paper_v4/results", showWarnings = FALSE, recursive = TRUE)

# Define colour schemes
age_colors <- c(Young = "#4575B4", Old = "#D73027")
age_colors_lc <- c(young = "#4575B4", old = "#D73027")
time_colors <- c("2" = "#FDE725", "10" = "#5DC863", "20" = "#21908C", "30" = "#440154")

# ================================================================
# LOAD DATA
# ================================================================

cat("\n============================================================\n")
cat("PART 1: LOADING DATA\n")
cat("============================================================\n")

# Load master files
master_counts <- read.csv("data/PRT4_paper_data/master_gene_counts.csv")
master_meta <- read.csv("data/PRT4_paper_data/master_metadata.csv")

cat("Master counts:", nrow(master_counts), "genes x", ncol(master_counts) - 2, "samples\n")
cat("Master metadata:", nrow(master_meta), "samples\n")

# Convert to matrix format
gene_symbols <- master_counts$gene_symbol
ensembl_ids <- master_counts$ensembl_gene_id
count_matrix <- as.matrix(master_counts[, -c(1, 2)])
rownames(count_matrix) <- gene_symbols

# Handle duplicate gene symbols
dup_genes <- duplicated(rownames(count_matrix))
if (sum(dup_genes) > 0) {
  cat("Removing", sum(dup_genes), "duplicate gene symbols\n")
  count_matrix <- count_matrix[!dup_genes, ]
}

cat("Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# Setup metadata
rownames(master_meta) <- master_meta$sample_id
stopifnot(all(colnames(count_matrix) %in% master_meta$sample_id))
master_meta <- master_meta[colnames(count_matrix), ]

cat("\nDataset breakdown:\n")
print(table(master_meta$dataset))
cat("\nAge group breakdown:\n")
print(table(master_meta$dataset, master_meta$age_group))
cat("\nTimepoint breakdown by dataset:\n")
print(table(master_meta$dataset, master_meta$timepoint, useNA = "ifany"))

# ================================================================
# SPLIT INTO ANALYSIS SUBSETS
# ================================================================

# PRT2 (Young)
meta_young <- master_meta %>%
  filter(dataset == "PRT2", muscle == "TA", stim == "Spillover") %>%
  mutate(
    time = factor(timepoint, levels = c(2, 10, 20, 30)),
    limb = factor(ifelse(limb == "Left", "Stim", "Control"), levels = c("Control", "Stim")),
    animal_id = factor(animal_id),
    sex = factor(sex)
  )

# PRT4 (Old)
meta_old <- master_meta %>%
  filter(dataset == "PRT4", muscle == "TA", stim == "Spillover", timepoint %in% c(2, 10, 20, 30)) %>%
  mutate(
    time = factor(timepoint, levels = c(2, 10, 20, 30)),
    limb = factor(ifelse(limb == "Left", "Stim", "Control"), levels = c("Control", "Stim")),
    animal_id = factor(animal_id),
    sex = factor(sex)
  )

# Subset count matrices
counts_young <- count_matrix[, meta_young$sample_id]
counts_old <- count_matrix[, meta_old$sample_id]

cat("\n=== ANALYSIS SUBSETS (before filtering) ===\n")
cat("Young (PRT2):", ncol(counts_young), "samples,", length(unique(meta_young$animal_id)), "animals\n")
cat("Old (PRT4):", ncol(counts_old), "samples,", length(unique(meta_old$animal_id)), "animals\n")

# Fix 432 column swap
cols_432 <- grep("432", colnames(counts_old), value = TRUE)
if (length(cols_432) == 2) {
  col_L <- cols_432[grepl("TAL", cols_432)]
  col_R <- cols_432[grepl("TAR", cols_432)]
  if (length(col_L) > 0 & length(col_R) > 0) {
    tmp <- counts_old[, col_L]
    counts_old[, col_L] <- counts_old[, col_R]
    counts_old[, col_R] <- tmp
    cat("Fixed 432 L/R swap\n")
  }
}

# ================================================================
# GENE FILTERING
# ================================================================

cat("\n============================================================\n")
cat("GENE FILTERING\n")
cat("============================================================\n")

min_count <- 10
min_samples_young <- ncol(counts_young) / 2
min_samples_old <- ncol(counts_old) / 2

keep_young <- rowSums(counts_young >= min_count) >= min_samples_young
keep_old <- rowSums(counts_old >= min_count) >= min_samples_old
keep_genes <- keep_young & keep_old

cat("Genes before filtering:", nrow(counts_young), "\n")
cat("Genes passing filter in Young:", sum(keep_young), "\n")
cat("Genes passing filter in Old:", sum(keep_old), "\n")
cat("Genes passing filter in BOTH:", sum(keep_genes), "\n")

counts_young <- counts_young[keep_genes, ]
counts_old <- counts_old[keep_genes, ]

cat("\n=== ANALYSIS SUBSETS (after filtering) ===\n")
cat("Young:", nrow(counts_young), "genes x", ncol(counts_young), "samples\n")
cat("Old:", nrow(counts_old), "genes x", ncol(counts_old), "samples\n")

# ================================================================
# COMBINED DATA
# ================================================================

meta_combined <- rbind(
  meta_young %>% mutate(Age_group = "Young"),
  meta_old %>% mutate(Age_group = "Old")
) %>%
  mutate(Age_group = factor(Age_group, levels = c("Young", "Old")))

counts_combined <- cbind(counts_young, counts_old)

cat("\nCombined dataset:", ncol(counts_combined), "samples\n")
cat("  Young:", sum(meta_combined$Age_group == "Young"), "\n")
cat("  Old:", sum(meta_combined$Age_group == "Old"), "\n")

# ================================================================
# SHAVLAKADZE DATA
# ================================================================

meta_shav <- master_meta %>% filter(dataset == "Shavlakadze2019")
shav_sra_path <- "data/PRT4_paper_data/Shavlakadze2019/SraRunTable.csv"

if (file.exists(shav_sra_path) && file.size(shav_sra_path) > 0) {
  shav_sra <- read.csv(shav_sra_path)
  cat("\nShavlakadze SRA metadata loaded:", nrow(shav_sra), "samples\n")
  
  run_col <- colnames(shav_sra)[grep("Run|run|SRR", colnames(shav_sra), ignore.case = TRUE)[1]]
  shav_age_info <- data.frame(sample_id = shav_sra[[run_col]], age_months = shav_sra$AGE)
  
  meta_shav <- meta_shav %>%
    left_join(shav_age_info, by = "sample_id") %>%
    mutate(
      age_group = case_when(
        age_months %in% c("6Mo", "9Mo") ~ "young",
        age_months %in% c("12Mo") ~ "middle",
        age_months %in% c("18Mo", "21Mo", "24Mo", "27Mo") ~ "old",
        TRUE ~ NA_character_
      )
    )
  
  cat("Age groups assigned:\n")
  print(table(meta_shav$age_group, meta_shav$age_months, useNA = "ifany"))
}

counts_shav <- count_matrix[keep_genes, meta_shav$sample_id]
cat("\nShavlakadze samples:", ncol(counts_shav), "\n")

# Legacy compatibility
sample_metadata <- meta_combined %>%
  dplyr::select(sample_id, animal_id, limb, sex, stim, Age_group, time) %>%
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$sample_id
sample_metadata$Age_group <- tolower(as.character(sample_metadata$Age_group))
merged_data <- counts_combined

cat("\n============================================================\n")
cat("DATA LOADING COMPLETE\n")
cat("============================================================\n")

# ============================================================
# FIGURE 1: PHENOTYPIC DATA - MUSCLE HYPERTROPHY
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 1: PHENOTYPIC DATA - MUSCLE HYPERTROPHY\n")
cat("============================================================\n")

# Load collection data
collection_path <- "data/PRT4_paper_data/collection_spreadsheet_131125.xlsx"
collection_data <- read_excel(collection_path)

# Get animal IDs with RNA-seq data
young_animals_with_rna <- as.character(unique(meta_young$animal_id))
old_animals_with_rna <- as.character(unique(meta_old$animal_id))

# Process hypertrophy data
hypertrophy_data <- collection_data %>%
  filter(stim == "Spillover") %>%
  filter(time %in% c("2", "10", "20", "30")) %>%
  mutate(
    animal_id = as.character(JFTO),
    TAL = as.numeric(TALwhole),
    TAR = as.numeric(TARwhole),
    body_weight = as.numeric(`harvestweight(g)`),
    TAL_norm = TAL / body_weight,
    TAR_norm = TAR / body_weight,
    pct_hypertrophy = (TAL - TAR) / TAR * 100,
    Time = factor(time, levels = c("2", "10", "20", "30")),
    Age = factor(ifelse(Age_group == "young", "Young", "Old"), levels = c("Young", "Old"))
  ) %>%
  filter(!is.na(TAL) & !is.na(TAR) & !is.na(body_weight)) %>%
  filter(
    (Age == "Young" & animal_id %in% young_animals_with_rna) |
      (Age == "Old" & animal_id %in% old_animals_with_rna)
  ) %>%
  dplyr::select(animal_id, Time, Age, TAL, TAR, body_weight, TAL_norm, TAR_norm, pct_hypertrophy)

cat("\nSample sizes:\n")
print(table(hypertrophy_data$Age, hypertrophy_data$Time))

# Reshape to long format
mass_long <- hypertrophy_data %>%
  pivot_longer(cols = c(TAL_norm, TAR_norm), names_to = "Stimulation", values_to = "mass_norm") %>%
  mutate(Stimulation = factor(ifelse(Stimulation == "TAL_norm", "Stimulated", "Control"), 
                              levels = c("Control", "Stimulated")))

# ============================================================
# STATISTICAL ANALYSIS
# ============================================================

# Panel A: LMM for muscle mass
lmm_mass <- lmer(mass_norm ~ Age * Stimulation * Time + (1|animal_id), data = mass_long)
mass_anova <- anova(lmm_mass, type = 3)

mass_p <- data.frame(
  Effect = rownames(mass_anova),
  F_value = round(mass_anova$`F value`, 2),
  p_value = mass_anova$`Pr(>F)`
)

cat("\nPanel A ANOVA:\n")
print(mass_p)

# Post-hoc for Age x Stimulation
cat("\nPost-hoc: Stimulation effect by Age:\n")
emm_stim <- emmeans(lmm_mass, pairwise ~ Stimulation | Age)
print(emm_stim$contrasts)

# Panel B: Linear model for hypertrophy
lm_pct <- lm(pct_hypertrophy ~ Age * Time, data = hypertrophy_data)
pct_anova <- anova(lm_pct)
age_p <- pct_anova$`Pr(>F)`[1]

cat("\nHypertrophy ANOVA:\n")
print(pct_anova)

# Summary statistics
mass_summary <- mass_long %>%
  group_by(Age, Time, Stimulation) %>%
  summarise(mean = mean(mass_norm), sd = sd(mass_norm), n = n(), .groups = "drop")

pct_summary <- hypertrophy_data %>%
  group_by(Age, Time) %>%
  summarise(
    mean = mean(pct_hypertrophy), 
    sd = sd(pct_hypertrophy), 
    se = sd(pct_hypertrophy) / sqrt(n()),
    n = n(), 
    .groups = "drop"
  ) %>%
  mutate(time_numeric = as.numeric(as.character(Time)))

# Day 30 summary
young_d30 <- pct_summary$mean[pct_summary$Age == "Young" & pct_summary$Time == "30"]
old_d30 <- pct_summary$mean[pct_summary$Age == "Old" & pct_summary$Time == "30"]
diff_d30 <- round(young_d30 - old_d30, 0)

cat("\nDay 30 growth:\n")
cat("  Young:", round(young_d30, 1), "%\n")
cat("  Old:", round(old_d30, 1), "%\n")
cat("  Difference:", diff_d30, "percentage points\n")

# ============================================================
# FORMAT P-VALUES
# ============================================================

format_p <- function(p) {
  if (p < 0.001) return("<0.001")
  else return(paste0("=", format(round(p, 3), nsmall = 3)))
}

mass_stats <- paste0(
  "LMM: Age p", format_p(mass_p$p_value[mass_p$Effect == "Age"]),
  ", Stim p", format_p(mass_p$p_value[mass_p$Effect == "Stimulation"]),
  ", Age×Stim p", format_p(mass_p$p_value[mass_p$Effect == "Age:Stimulation"]),
  ", Stim×Time p", format_p(mass_p$p_value[mass_p$Effect == "Stimulation:Time"])
)

cat("\nPanel A annotation:", mass_stats, "\n")

# ============================================================
# PANEL A: MUSCLE MASS (FACETED)
# ============================================================

fig1a <- ggplot(mass_summary, aes(x = as.numeric(as.character(Time)), y = mean * 1000,
                                  color = Stimulation, group = Stimulation)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (mean - sd) * 1000, ymax = (mean + sd) * 1000), 
                width = 1.2, linewidth = 0.6) +
  geom_point(data = mass_long %>% mutate(mass_mg = mass_norm * 1000),
             aes(x = as.numeric(as.character(Time)), y = mass_mg, color = Stimulation),
             alpha = 0.4, size = 2,
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1)) +
  facet_wrap(~ Age) +
  scale_color_manual(values = c("Control" = "gray50", "Stimulated" = "steelblue"), name = "") +
  scale_x_continuous(breaks = c(2, 10, 20, 30)) +
  labs(x = "Day", y = "TA Mass / Body Weight (mg/g)", 
       title = "A: Muscle Mass",
       subtitle = mass_stats) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, face = "italic", color = "gray30")
  )

# ============================================================
# PANEL B: HYPERTROPHY TRAJECTORY
# ============================================================

fig1b <- ggplot(pct_summary, aes(x = time_numeric, y = mean, color = Age, group = Age)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1, linewidth = 0.8) +
  geom_point(data = hypertrophy_data %>% mutate(time_numeric = as.numeric(as.character(Time))),
             aes(x = time_numeric, y = pct_hypertrophy, color = Age),
             alpha = 0.4, size = 2, position = position_jitter(width = 0.5)) +
  scale_color_manual(values = age_colors) +
  scale_x_continuous(breaks = c(2, 10, 20, 30)) +
  labs(
    x = "Day", 
    y = "Muscle Growth (%)",
    title = "B: Hypertrophy Over Time",
    subtitle = paste0("Age p = ", signif(age_p, 2), 
                      "; Day 30: Young ", round(young_d30, 0), "%, Old ", round(old_d30, 0), "%")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic", color = "gray30")
  )

# ============================================================
# COMBINE AND SAVE FIGURE 1
# ============================================================

fig1 <- (fig1a | fig1b) + plot_layout(widths = c(1.5, 1))

ggsave("paper_v4/figures/Figure1_hypertrophy.png", fig1, width = 14, height = 5, dpi = 300, bg = "white")

cat("\nFigure 1 saved.\n")

# Save results
write.csv(hypertrophy_data, "paper_v4/results/hypertrophy_individual_data.csv", row.names = FALSE)
write.csv(pct_summary, "paper_v4/results/hypertrophy_summary.csv", row.names = FALSE)
write.csv(mass_p, "paper_v4/results/lmm_mass_anova.csv", row.names = FALSE)

cat("\n=== FIGURE 1 COMPLETE ===\n")

# ================================================================
# FIGURE 2: TRANSCRIPTOME OVERVIEW
# ================================================================

cat("\n============================================================\n")
cat("FIGURE 2: TRANSCRIPTOME OVERVIEW\n")
cat("============================================================\n")

# ============================================================
# PCA ANALYSIS
# ============================================================

run_pca <- function(counts, meta) {
  dge <- calcNormFactors(DGEList(counts = counts))
  logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  pca <- prcomp(t(logcpm), scale. = TRUE)
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  scores <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta, by = "sample_id")
  list(scores = scores, var = var_exp)
}

pca_young <- run_pca(counts_young, meta_young)
pca_old <- run_pca(counts_old, meta_old)

# Flip PC1 for consistent orientation
if (mean(pca_young$scores$PC1[pca_young$scores$limb == "Control"]) < 0) {
  pca_young$scores$PC1 <- -pca_young$scores$PC1
}

ctrl_mean <- mean(pca_old$scores$PC1[pca_old$scores$limb == "Control"])
stim_mean <- mean(pca_old$scores$PC1[pca_old$scores$limb == "Stim"])
if (stim_mean > ctrl_mean) {
  pca_old$scores$PC1 <- -pca_old$scores$PC1
}

# Set consistent axis limits
pc1_range <- range(c(pca_young$scores$PC1, pca_old$scores$PC1))
pc2_range <- range(c(pca_young$scores$PC2, pca_old$scores$PC2))
pc1_lim <- c(floor(pc1_range[1]/25)*25 - 25, ceiling(pc1_range[2]/25)*25 + 25)
pc2_lim <- c(floor(pc2_range[1]/25)*25 - 25, ceiling(pc2_range[2]/25)*25 + 25)

# Panel A: Young PCA
fig2a <- ggplot(pca_young$scores, aes(PC1, PC2, color = time, shape = limb)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = time_colors) +
  scale_shape_manual(values = c(Control = 16, Stim = 17)) +
  coord_cartesian(xlim = pc1_lim, ylim = pc2_lim) +
  labs(x = paste0("PC1 (", pca_young$var[1], "%)"), 
       y = paste0("PC2 (", pca_young$var[2], "%)"), 
       title = "A: Young") +
  theme_classic(base_size = 11) + 
  theme(legend.position = "right")

# Panel B: Old PCA
fig2b <- ggplot(pca_old$scores, aes(PC1, PC2, color = time, shape = limb)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = time_colors) +
  scale_shape_manual(values = c(Control = 16, Stim = 17)) +
  coord_cartesian(xlim = pc1_lim, ylim = pc2_lim) +
  labs(x = paste0("PC1 (", pca_old$var[1], "%)"), 
       y = paste0("PC2 (", pca_old$var[2], "%)"), 
       title = "B: Old") +
  theme_classic(base_size = 11) + 
  theme(legend.position = "right")

# ============================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================

run_de <- function(counts, meta, include_sex = FALSE) {
  if (include_sex) {
    design <- model.matrix(~ sex + time * limb, data = meta)
  } else {
    design <- model.matrix(~ time * limb, data = meta)
  }
  colnames(design) <- make.names(colnames(design))
  
  dge <- calcNormFactors(DGEList(counts = counts))
  v <- voom(dge, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = meta$animal_id)
  fit <- lmFit(v, design, block = meta$animal_id, correlation = corfit$consensus)
  fit <- eBayes(fit, robust = TRUE)
  
  contrasts <- makeContrasts(
    Day2 = limbStim, 
    Day10 = limbStim + time10.limbStim,
    Day20 = limbStim + time20.limbStim, 
    Day30 = limbStim + time30.limbStim,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, contrasts), robust = TRUE)
  
  cat("  Within-animal correlation:", round(corfit$consensus, 3), "\n")
  
  list(
    Day2 = topTable(fit2, coef = "Day2", number = Inf),
    Day10 = topTable(fit2, coef = "Day10", number = Inf),
    Day20 = topTable(fit2, coef = "Day20", number = Inf),
    Day30 = topTable(fit2, coef = "Day30", number = Inf)
  )
}

cat("\nYoung:\n")
results_young <- run_de(counts_young, meta_young, include_sex = FALSE)

cat("\nOld (with sex covariate):\n")
results_old <- run_de(counts_old, meta_old, include_sex = TRUE)

# ============================================================
# DEG COUNTS
# ============================================================

deg_counts <- data.frame(
  Day = rep(c("2", "10", "20", "30"), 2),
  Age = rep(c("Young", "Old"), each = 4),
  Up = c(
    sum(results_young$Day2$adj.P.Val < 0.05 & results_young$Day2$logFC > 0),
    sum(results_young$Day10$adj.P.Val < 0.05 & results_young$Day10$logFC > 0),
    sum(results_young$Day20$adj.P.Val < 0.05 & results_young$Day20$logFC > 0),
    sum(results_young$Day30$adj.P.Val < 0.05 & results_young$Day30$logFC > 0),
    sum(results_old$Day2$adj.P.Val < 0.05 & results_old$Day2$logFC > 0),
    sum(results_old$Day10$adj.P.Val < 0.05 & results_old$Day10$logFC > 0),
    sum(results_old$Day20$adj.P.Val < 0.05 & results_old$Day20$logFC > 0),
    sum(results_old$Day30$adj.P.Val < 0.05 & results_old$Day30$logFC > 0)
  ),
  Down = c(
    sum(results_young$Day2$adj.P.Val < 0.05 & results_young$Day2$logFC < 0),
    sum(results_young$Day10$adj.P.Val < 0.05 & results_young$Day10$logFC < 0),
    sum(results_young$Day20$adj.P.Val < 0.05 & results_young$Day20$logFC < 0),
    sum(results_young$Day30$adj.P.Val < 0.05 & results_young$Day30$logFC < 0),
    sum(results_old$Day2$adj.P.Val < 0.05 & results_old$Day2$logFC < 0),
    sum(results_old$Day10$adj.P.Val < 0.05 & results_old$Day10$logFC < 0),
    sum(results_old$Day20$adj.P.Val < 0.05 & results_old$Day20$logFC < 0),
    sum(results_old$Day30$adj.P.Val < 0.05 & results_old$Day30$logFC < 0)
  )
) %>%
  mutate(
    Day = factor(Day, levels = c("2", "10", "20", "30")),
    Age = factor(Age, levels = c("Young", "Old"))
  )

cat("\n=== DEG SUMMARY ===\n")
print(deg_counts)

# Panel C: DEG bar plot
deg_long <- deg_counts %>%
  pivot_longer(c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Count = ifelse(Direction == "Down", -Count, Count),
    Direction = factor(Direction, levels = c("Up", "Down"))
  )

fig2c <- ggplot(deg_long, aes(x = Day, y = Count, fill = Direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_wrap(~Age) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  scale_y_continuous(labels = abs) +
  labs(x = "Day", y = "Number of DEGs (FDR < 0.05)", title = "C: Differentially Expressed Genes") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey90"))

# ============================================================
# QUADRANT PLOTS (YOUNG VS OLD CORRELATION)
# ============================================================

make_quadrant <- function(y_res, o_res, day_label) {
  df <- merge(
    data.frame(gene = rownames(y_res), lfc_y = y_res$logFC, p_y = y_res$adj.P.Val),
    data.frame(gene = rownames(o_res), lfc_o = o_res$logFC, p_o = o_res$adj.P.Val),
    by = "gene"
  ) %>%
    mutate(
      cat = case_when(
        p_y < 0.05 & p_o < 0.05 ~ "Both",
        p_y < 0.05 ~ "Young only",
        p_o < 0.05 ~ "Old only",
        TRUE ~ "NS"
      ),
      cat = factor(cat, levels = c("Both", "Young only", "Old only", "NS"))
    )
  
  r <- cor(df$lfc_y, df$lfc_o, use = "complete.obs")
  
  ggplot(df, aes(x = lfc_y, y = lfc_o, color = cat)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, linetype = "dotted", color = "gray30") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(
      values = c("Both" = "orange", "Young only" = "#4575B4", "Old only" = "#D73027", "NS" = "gray80"),
      name = "DEG Status"
    ) +
    coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
    annotate("text", x = -5, y = 5, label = paste0("r = ", round(r, 2)), fontface = "plain") +
    labs(title = day_label, x = "Young log2FC", y = "Old log2FC") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
}

q2 <- make_quadrant(results_young$Day2, results_old$Day2, "Day 2")
q10 <- make_quadrant(results_young$Day10, results_old$Day10, "Day 10")
q20 <- make_quadrant(results_young$Day20, results_old$Day20, "Day 20")
q30 <- make_quadrant(results_young$Day30, results_old$Day30, "Day 30")

# Create legend
legend_df <- data.frame(
  x = 1:4, y = 1:4,
  cat = factor(c("Both", "Young only", "Old only", "NS"), 
               levels = c("Both", "Young only", "Old only", "NS"))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = cat)) +
  geom_point() +
  scale_color_manual(
    values = c("Both" = "orange", "Young only" = "#4575B4", "Old only" = "#D73027", "NS" = "gray80"),
    name = "DEG Status"
  ) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 9))

legend_grob <- cowplot::get_legend(legend_plot)

fig2d_plots <- (q2 | q10 | q20 | q30) / wrap_elements(legend_grob) + plot_layout(heights = c(10, 1))
fig2d <- wrap_elements(full = fig2d_plots) + 
  ggtitle("D: Transcriptional Correlation (Young vs Old log2FC)") +
  theme(plot.title = element_text(size = 11))

# ============================================================
# COMBINE AND SAVE FIGURE 2
# ============================================================

fig2 <- (fig2a | fig2b) / fig2c / fig2d + plot_layout(heights = c(1, 0.8, 1))

ggsave("paper_v4/figures/Figure2_transcriptome_overview.png", fig2, width = 14, height = 16, dpi = 300)

cat("Figure 2 saved\n")

# Save results
write.csv(deg_counts, "paper_v4/results/deg_counts_summary.csv", row.names = FALSE)

for (tp in c("Day2", "Day10", "Day20", "Day30")) {
  write.csv(results_young[[tp]] %>% rownames_to_column("gene"), 
            paste0("paper_v4/results/DE_young_", tp, ".csv"), row.names = FALSE)
  write.csv(results_old[[tp]] %>% rownames_to_column("gene"), 
            paste0("paper_v4/results/DE_old_", tp, ".csv"), row.names = FALSE)
}

cat("Individual DE results saved\n")

# ============================================================
# DEG OVERLAP ANALYSIS (VENN DIAGRAMS)
# ============================================================

cat("\n=== DEG OVERLAP ANALYSIS ===\n")

deg_young_d2 <- rownames(results_young$Day2)[results_young$Day2$adj.P.Val < 0.05]
deg_young_d10 <- rownames(results_young$Day10)[results_young$Day10$adj.P.Val < 0.05]
deg_young_d20 <- rownames(results_young$Day20)[results_young$Day20$adj.P.Val < 0.05]
deg_young_d30 <- rownames(results_young$Day30)[results_young$Day30$adj.P.Val < 0.05]

deg_old_d2 <- rownames(results_old$Day2)[results_old$Day2$adj.P.Val < 0.05]
deg_old_d10 <- rownames(results_old$Day10)[results_old$Day10$adj.P.Val < 0.05]
deg_old_d20 <- rownames(results_old$Day20)[results_old$Day20$adj.P.Val < 0.05]
deg_old_d30 <- rownames(results_old$Day30)[results_old$Day30$adj.P.Val < 0.05]

young_list <- list(D2 = deg_young_d2, D10 = deg_young_d10, D20 = deg_young_d20, D30 = deg_young_d30)
old_list <- list(D2 = deg_old_d2, D10 = deg_old_d10, D20 = deg_old_d20, D30 = deg_old_d30)

venn_young_gg <- ggVennDiagram(young_list, label_alpha = 0, 
                               set_color = c("#FDE725", "#5DC863", "#21908C", "#440154")) +
  scale_fill_gradient(low = "white", high = "#4575B4") +
  labs(title = "Young DEGs") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

venn_old_gg <- ggVennDiagram(old_list, label_alpha = 0,
                             set_color = c("#FDE725", "#5DC863", "#21908C", "#440154")) +
  scale_fill_gradient(low = "white", high = "#D73027") +
  labs(title = "Old DEGs") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

venn_combined <- venn_young_gg | venn_old_gg
ggsave("paper_v4/figures/FigureS1_DEG_overlap_venn.png", venn_combined, width = 16, height = 8, dpi = 150)

# Overlap statistics
young_all_tp <- Reduce(intersect, list(deg_young_d2, deg_young_d10, deg_young_d20, deg_young_d30))
old_all_tp <- Reduce(intersect, list(deg_old_d2, deg_old_d10, deg_old_d20, deg_old_d30))
young_any_tp <- Reduce(union, list(deg_young_d2, deg_young_d10, deg_young_d20, deg_young_d30))
old_any_tp <- Reduce(union, list(deg_old_d2, deg_old_d10, deg_old_d20, deg_old_d30))

cat("\nYoung DEG overlap:\n")
cat("  DEGs in all 4 timepoints:", length(young_all_tp), "\n")
cat("  DEGs in any timepoint:", length(young_any_tp), "\n")

cat("  Proportion persistent:", round(length(young_all_tp)/length(young_any_tp)*100, 1), "%\n")

cat("\nOld DEG overlap:\n")
cat("  DEGs in all 4 timepoints:", length(old_all_tp), "\n")
cat("  DEGs in any timepoint:", length(old_any_tp), "\n")
cat("  Proportion persistent:", round(length(old_all_tp)/length(old_any_tp)*100, 1), "%\n")

# Save overlap stats
deg_overlap_stats <- data.frame(
  Age = c("Young", "Young", "Young", "Young", "Old", "Old", "Old", "Old"),
  Metric = rep(c("D2", "D10", "D20", "D30"), 2),
  n_DEGs = c(length(deg_young_d2), length(deg_young_d10), length(deg_young_d20), length(deg_young_d30),
             length(deg_old_d2), length(deg_old_d10), length(deg_old_d20), length(deg_old_d30))
)

write.csv(deg_overlap_stats, "paper_v4/results/deg_overlap_stats.csv", row.names = FALSE)
write.csv(data.frame(gene = young_all_tp), "paper_v4/results/core_degs_young.csv", row.names = FALSE)
write.csv(data.frame(gene = old_all_tp), "paper_v4/results/core_degs_old.csv", row.names = FALSE)

cat("\n=== FIGURE 2 COMPLETE ===\n")

# ============================================================
# FIGURE 3: GSEA PATHWAY DIVERGENCE
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 3: GSEA PATHWAY DIVERGENCE\n")
cat("============================================================\n")

# Get Hallmark gene sets
hallmark_sets <- msigdbr(species = "Rattus norvegicus", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

cat("Loaded", length(hallmark_list), "Hallmark pathways\n")

# ============================================================
# RUN GSEA FOR EACH AGE x TIMEPOINT
# ============================================================

gsea_results_list <- list()

for (age in c("Young", "Old")) {
  results_list <- if (age == "Young") results_young else results_old
  
  for (day in c("Day2", "Day10", "Day20", "Day30")) {
    res <- results_list[[day]]
    
    # Create ranked list
    ranked <- res$t
    names(ranked) <- rownames(res)
    ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)
    
    # Run GSEA
    set.seed(42)
    gsea_res <- fgsea(hallmark_list, ranked, minSize = 15, maxSize = 500, nPermSimple = 10000)
    
    gsea_res$Age <- age
    gsea_res$Day <- day
    gsea_results_list[[paste(age, day, sep = "_")]] <- gsea_res
    
    cat("  ", age, day, ":", sum(gsea_res$padj < 0.05), "significant pathways\n")
  }
}

# Combine results
gsea_all <- bind_rows(gsea_results_list)

# Clean pathway names
gsea_all <- gsea_all %>%
  mutate(
    pathway = gsub("HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway),
    Day_num = as.numeric(gsub("Day", "", Day))
  )

# ============================================================
# CREATE WIDE FORMAT FOR COMPARISON
# ============================================================

hallmark_combined <- gsea_all %>%
  dplyr::select(pathway, Day, Age, NES, padj) %>%
  pivot_wider(
    names_from = Age,
    values_from = c(NES, padj),
    names_sep = "_"
  ) %>%
  rename(
    NES_young = NES_Young,
    NES_old = NES_Old,
    padj_young = padj_Young,
    padj_old = padj_Old
  ) %>%
  mutate(
    NES_diff = NES_old - NES_young,
    Day_num = as.numeric(gsub("Day", "", Day))
  )

# ============================================================
# SELECT TOP 15 DIVERGENT PATHWAYS
# ============================================================

pathway_divergence <- hallmark_combined %>%
  group_by(pathway) %>%
  summarise(
    max_abs_NES_diff = max(abs(NES_diff), na.rm = TRUE),
    mean_NES_diff = mean(NES_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_NES_diff))

top15_divergent <- pathway_divergence$pathway[1:15]

# ============================================================
# PREPARE DATA FOR PLOTTING
# ============================================================

plot_data <- gsea_all %>%
  filter(pathway %in% top15_divergent) %>%
  mutate(
    pathway = factor(pathway, levels = top15_divergent),
    significant = padj < 0.05
  )

# ============================================================
# FIGURE 3: TRAJECTORY PLOTS (3x5 GRID)
# ============================================================

fig3 <- ggplot(plot_data, aes(x = Day_num, y = NES, color = Age, group = Age)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = significant, fill = Age), size = 3, stroke = 1) +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21), guide = "none") +
  geom_point(
    data = plot_data %>% filter(!significant),
    aes(x = Day_num, y = NES),
    shape = 21, size = 3, stroke = 1, fill = "white"
  ) +
  facet_wrap(~ pathway, ncol = 5, scales = "free_y") +
  scale_x_continuous(breaks = c(2, 10, 20, 30), labels = c("2", "10", "20", "30")) +
  labs(
    x = "Day",
    y = "NES",
    title = "Top 15 Divergent Hallmark Pathways",
    subtitle = "Filled = FDR < 0.05, Open = FDR >= 0.05"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  )

# Save
ggsave("paper_v4/figures/Figure3_GSEA.png", fig3, width = 14, height = 10, dpi = 300, bg = "white")

cat("\n=== FIGURE 3 COMPLETE ===\n")

# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(gsea_all %>% dplyr::select(-leadingEdge), 
          "paper_v4/results/gsea_all_results.csv", row.names = FALSE)
write.csv(pathway_divergence, 
          "paper_v4/results/pathway_divergence_metrics.csv", row.names = FALSE)
write.csv(hallmark_combined, 
          "paper_v4/results/hallmark_combined_young_old.csv", row.names = FALSE)

cat("GSEA results saved\n")

# ================================================================
# FIGURE 4: sPLS-DA WITH ENRICHMENT
# ================================================================

cat("\n============================================================\n")
cat("FIGURE 4: sPLS-DA WITH ENRICHMENT\n")
cat("============================================================\n")

# ============================================================
# CREATE LOG2FC MATRIX
# ============================================================

meta_all <- meta_combined %>%
  mutate(
    time = factor(time, levels = c("2", "10", "20", "30")),
    limb = factor(limb, levels = c("Control", "Stim"))
  )

counts_all <- counts_combined[, meta_all$sample_id]
dge_all <- calcNormFactors(DGEList(counts = counts_all))
v_all <- voom(dge_all, model.matrix(~ time * limb, data = meta_all), plot = FALSE)

# Calculate log2FC for each animal (Stim - Control)
log2FC_list <- list()
for (aid in unique(meta_all$animal_id)) {
  sub <- meta_all %>% filter(animal_id == aid)
  if (nrow(sub) != 2) next
  s <- sub$sample_id[sub$limb == "Stim"]
  c <- sub$sample_id[sub$limb == "Control"]
  if (all(c(s, c) %in% colnames(v_all$E))) {
    log2FC_list[[as.character(aid)]] <- v_all$E[, s] - v_all$E[, c]
  }
}

log2FC_mat <- do.call(cbind, log2FC_list)
log2FC_mat[!is.finite(log2FC_mat)] <- 0

cat("log2FC matrix:", nrow(log2FC_mat), "genes x", ncol(log2FC_mat), "animals\n")

# ============================================================
# PREPARE METADATA FOR sPLS-DA
# ============================================================

mint_meta <- meta_all %>%
  dplyr::select(animal_id, Age_group, time, sex) %>%
  distinct() %>%
  mutate(animal_id = as.character(animal_id)) %>%
  filter(animal_id %in% colnames(log2FC_mat))

mint_meta <- mint_meta[match(colnames(log2FC_mat), mint_meta$animal_id), ]

cat("Animals in analysis:", nrow(mint_meta), "\n")
cat("Age distribution:\n")
print(table(mint_meta$Age_group))

# ============================================================
# RUN sPLS-DA
# ============================================================

X <- t(log2FC_mat)
Y <- factor(mint_meta$Age_group, levels = c("Young", "Old"))

cat("\nY levels:", levels(Y), "\n")
cat("Y table:\n")
print(table(Y))

stopifnot(nrow(X) == length(Y))
stopifnot(all(!is.na(Y)))

set.seed(1000)
cat("\nRunning sPLS-DA with keepX = c(100, 200)...\n")

splsda_fit <- splsda(X = X, Y = Y, ncomp = 2, keepX = c(200))

cat("sPLS-DA complete\n")
cat("Component 1 variance:", round(splsda_fit$prop_expl_var$X[1] * 100, 1), "%\n")
cat("Component 2 variance:", round(splsda_fit$prop_expl_var$X[2] * 100, 1), "%\n")

# ============================================================
# EXTRACT SCORES AND LOADINGS
# ============================================================

df_scores <- data.frame(
  animal_id = mint_meta$animal_id,
  Comp1 = splsda_fit$variates$X[, 1],
  Comp2 = splsda_fit$variates$X[, 2],
  Age = Y,
  Time = factor(mint_meta$time)
)

# Determine component direction
young_c1_mean <- mean(df_scores$Comp1[df_scores$Age == "Young"])
old_c1_mean <- mean(df_scores$Comp1[df_scores$Age == "Old"])

if (young_c1_mean > old_c1_mean) {
  pos_direction <- "Young"
  neg_direction <- "Old"
} else {
  pos_direction <- "Old"
  neg_direction <- "Young"
}

cat("Component 1 direction: Positive =", pos_direction, ", Negative =", neg_direction, "\n")

# Extract loadings
comp1_loadings <- data.frame(
  gene = rownames(splsda_fit$loadings$X),
  loading = as.numeric(splsda_fit$loadings$X[, 1])
) %>%
  filter(loading != 0) %>%
  arrange(desc(loading)) %>%
  mutate(direction = ifelse(loading > 0, pos_direction, neg_direction))

cat("Total genes with non-zero loadings:", nrow(comp1_loadings), "\n")
cat("Young-associated:", sum(comp1_loadings$direction == "Young"), "\n")
cat("Old-associated:", sum(comp1_loadings$direction == "Old"), "\n")

# ============================================================
# PANEL A: SCORE PLOT
# ============================================================

fig4a <- ggplot(df_scores, aes(x = Comp1, y = Comp2, color = Age, shape = Time)) +
  geom_point(size = 3.5, alpha = 0.8) +
  stat_ellipse(aes(group = Age), level = 0.95, linewidth = 1) +
  scale_color_manual(values = age_colors) +
  scale_shape_manual(values = c("2" = 16, "10" = 17, "20" = 15, "30" = 3)) +
  labs(
    x = paste0("Component 1 (", round(splsda_fit$prop_expl_var$X[1] * 100, 1), "%)"),
    y = paste0("Component 2 (", round(splsda_fit$prop_expl_var$X[2] * 100, 1), "%)"),
    title = "A: sPLS-DA Score Plot"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

# ============================================================
# PANEL B: COMPONENT 1 TRAJECTORY
# ============================================================

fig4b <- ggplot(df_scores, aes(x = as.numeric(as.character(Time)), y = Comp1, color = Age, group = Age)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 1.5, linewidth = 0.8) +
  scale_color_manual(values = age_colors) +
  scale_x_continuous(breaks = c(2, 10, 20, 30)) +
  labs(x = "Day", y = "Component 1 Score", title = "B: Discrimination Over Time") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

# ============================================================
# PANEL C: TOP 25 LOADINGS
# ============================================================

top25 <- comp1_loadings %>%
  arrange(desc(abs(loading))) %>%
  head(25)

fig4c <- ggplot(top25, aes(x = loading, y = reorder(gene, loading), fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = c("Old" = "#D73027", "Young" = "#4575B4")) +
  labs(x = "Component 1 Loading", y = "", title = "C: Top 25 Discriminating Genes") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

# ============================================================
# PANEL D: GO:BP ENRICHMENT
# ============================================================

cat("\nRunning GO:BP enrichment...\n")

young_genes <- comp1_loadings %>% filter(direction == "Young") %>% pull(gene)
old_genes <- comp1_loadings %>% filter(direction == "Old") %>% pull(gene)

# Young enrichment
young_enrich <- enrichGO(
  gene = young_genes,
  OrgDb = org.Rn.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

n_young_terms <- ifelse(!is.null(young_enrich), nrow(as.data.frame(young_enrich)), 0)
cat("Young enrichment results:", n_young_terms, "terms\n")

# Old enrichment
old_enrich <- enrichGO(
  gene = old_genes,
  OrgDb = org.Rn.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

n_old_terms <- ifelse(!is.null(old_enrich), nrow(as.data.frame(old_enrich)), 0)
cat("Old enrichment results:", n_old_terms, "terms\n")

# Create enrichment plot
if (n_old_terms > 0) {
  cat("\nTop GO:BP terms (Old-associated genes):\n")
  print(head(as.data.frame(old_enrich)[, c("Description", "p.adjust", "Count")], 10))
  
  old_enrich_sim <- pairwise_termsim(old_enrich)
  
  fig4d <- emapplot(old_enrich_sim, showCategory = 20, color = "p.adjust", layout = "kk")
  
  # Remove text layers
  fig4d$layers <- fig4d$layers[!sapply(fig4d$layers, function(l) {
    "GeomText" %in% class(l$geom) || "GeomTextRepel" %in% class(l$geom) || "GeomLabel" %in% class(l$geom)
  })]
  
  fig4d <- fig4d +
    scale_fill_gradient(low = "#D73027", high = "#FDAE61", name = "p.adjust") +
    labs(
      title = "D: GO:BP Enrichment (Old-Associated)",
      subtitle = "Glycogen and lipid metabolism"
    ) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5, color = "gray30")
    )
  
  write.csv(as.data.frame(old_enrich), "paper_v4/results/splsda_old_GO_BP_full.csv", row.names = FALSE)
  
} else {
  fig4d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No significant GO:BP enrichment", size = 5) +
    theme_void() +
    ggtitle("D: GO:BP Enrichment")
}

# ============================================================
# COMBINE AND SAVE FIGURE 4
# ============================================================

fig4_top <- (fig4a | fig4b) + plot_layout(widths = c(1, 1))
fig4_bottom <- (fig4c | fig4d) + plot_layout(widths = c(0.8, 1))
fig4 <- fig4_top / fig4_bottom + plot_layout(heights = c(1, 1))

ggsave("paper_v4/figures/Figure4_sPLSDA.png", fig4, width = 14, height = 12, dpi = 300, bg = "white")

write.csv(comp1_loadings, "paper_v4/results/splsda_comp1_loadings_200.csv", row.names = FALSE)

cat("\n=== FIGURE 4 COMPLETE ===\n")

# ============================================================
# FIGURE 4 SUMMARY
# ============================================================

cat("\n=== FIGURE 4 SUMMARY ===\n")
cat("\nPanel A: sPLS-DA score plot\n")
cat("  Clear separation between Young and Old\n")
cat("  Component 1 explains", round(splsda_fit$prop_expl_var$X[1] * 100, 1), "% variance\n")

cat("\nPanel B: Discrimination over time\n")
cat("  Shows how age separation changes across timepoints\n")

cat("\nPanel C: Top 25 discriminating genes\n")
cat("  Young-associated (positive):", sum(top25$direction == "Young"), "genes\n")
cat("  Old-associated (negative):", sum(top25$direction == "Old"), "genes\n")

cat("\nPanel D: GO:BP enrichment\n")
if (n_old_terms > 0) {
  cat("  Old genes enriched for METABOLIC REGULATION\n")
} else {
  cat("  No significant enrichment\n")
}

cat("\n=== KEY FINDING ===\n")
cat("Old-associated genes converge on METABOLIC REGULATION\n")
cat("Young-associated genes show NO coherent functional signature\n")

# ============================================================
# FIGURE 5: MYC EXECUTION FAILURE
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 5: MYC EXECUTION FAILURE\n")
cat("============================================================\n")

# Get MYC target gene sets
hallmark_sets <- msigdbr(species = "Rattus norvegicus", category = "H")

myc_v1_genes <- hallmark_sets %>%
  filter(gs_name == "HALLMARK_MYC_TARGETS_V1") %>%
  pull(gene_symbol) %>%
  unique()

myc_v2_genes <- hallmark_sets %>%
  filter(gs_name == "HALLMARK_MYC_TARGETS_V2") %>%
  pull(gene_symbol) %>%
  unique()

cat("MYC V1 genes:", length(myc_v1_genes), "\n")
cat("MYC V2 genes:", length(myc_v2_genes), "\n")
cat("Overlap:", length(intersect(myc_v1_genes, myc_v2_genes)), "\n")

# ============================================================
# FUNCTIONAL CATEGORY ASSIGNMENT - SIMPLIFIED (4 CATEGORIES)
# ============================================================

assign_functional_category <- function(gene) {
  gene_upper <- toupper(gene)
  
  # Ribosome biogenesis (includes RNA processing, translation machinery, chaperones)
  if (grepl("^RPL[0-9]|^RPS[0-9]|^MRPL|^MRPS|^NOP[0-9]|^RRP|^DDX|^DHX|^LSM[0-9]|^SNRP|^UTP|^IMP[0-9]|^NHP2|^GAR1|^FBL$|^DKC1|^NOP56|^NOP58|^PPAN|^NOC|^EMG1|^WDR|^BOP1|^PES1|^EIF|^EEF|^PAIP|^RPP|^CCT[0-9]|^TCP1|^HSP[ABDE]|^HSPA|^HSPB|^HSPD|^HSPE|^HSPH|^DNAJ|^PPID|^PPIA$|^PPIB$|^PPIF|^FKBP|^CALR|^CANX|^PDIA|^ERP[0-9]", gene_upper)) {
    return("Ribosome biogenesis")
  }
  
  # Cell proliferation (cell cycle, DNA replication, nucleotide synthesis)
  if (grepl("^MCM[0-9]|^CDC[0-9]|^CDK[0-9]|^CCNA|^CCNB|^CCND|^CCNE|^E2F[0-9]|^ORC[0-9]|^MAD2|^BUB|^PCNA|^CDT1|^DBF4|^PLK|^AURK|^TOP2|^RFC[0-9]|^POLA|^POLD|^POLE|^PRIM|^GINS|^FEN1|^LIG1|^RRM[0-9]|^PRPS[0-9]|^IMPDH|^CTPS|^GMPS|^TYMS|^ODC1|^DCTPP|^CAD$|^UMPS|^ATIC|^GART|^PAICS|^PFAS|^ADSL|^ADSS|^PPAT|^NME[0-9]|^TK1$|^DTYMK|^UCK|^DCK$|^UNG", gene_upper)) {
    return("Cell proliferation")
  }
  
  # Metabolism
  if (grepl("^SLC[0-9]|^HK[0-9]|^LDHA|^LDHB|^ENO[0-9]|^PKM|^PGK[0-9]|^ALDOA|^ALDOB|^ALDOC|^TPI1|^GAPDH|^PFK|^GPI$|^PFKP|^PFKM|^PFKL|^PDH|^IDH|^OGDH|^SDHA|^SDHB|^MDH|^CS$|^ACO[0-9]|^FH$", gene_upper)) {
    return("Metabolism")
  }
  
  # Other (transcription, chromatin, signalling, everything else)
  return("Other")
}

# Define category colours
category_colors <- c(
  "Ribosome biogenesis" = "#1b9e77",
  "Cell proliferation" = "#d95f02",
  "Metabolism" = "#e6ab02",
  "Other" = "#666666"
)

# Categorise all MYC target genes
myc_gene_category <- data.frame(
  gene = unique(c(myc_v1_genes, myc_v2_genes))
) %>%
  mutate(
    in_v1 = gene %in% myc_v1_genes,
    in_v2 = gene %in% myc_v2_genes,
    functional_category = sapply(gene, assign_functional_category),
    functional_category = factor(functional_category, 
                                 levels = c("Ribosome biogenesis",
                                            "Cell proliferation",
                                            "Metabolism",
                                            "Other"))
  )

cat("\nFunctional categories:\n")
print(table(myc_gene_category$functional_category))

# ============================================================
# PANEL A: MYC EXPRESSION
# ============================================================

myc_gene <- "Myc"

# Get individual animal MYC log2FC
myc_individual <- data.frame()

for (aid in colnames(log2FC_mat)) {
  animal_meta <- mint_meta %>% filter(animal_id == aid)
  if (nrow(animal_meta) == 1 && myc_gene %in% rownames(log2FC_mat)) {
    myc_individual <- rbind(myc_individual, data.frame(
      animal_id = aid,
      Age = as.character(animal_meta$Age_group),
      Day = as.numeric(as.character(animal_meta$time)),
      MYC_logFC = log2FC_mat[myc_gene, aid]
    ))
  }
}

# Summary for plotting
myc_summary <- myc_individual %>%
  group_by(Age, Day) %>%
  summarise(
    mean_logFC = mean(MYC_logFC),
    se_logFC = sd(MYC_logFC) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Add significance from DE results
myc_sig <- data.frame()
for (age in c("Young", "Old")) {
  results_list <- if (age == "Young") results_young else results_old
  for (day in c("Day2", "Day10", "Day20", "Day30")) {
    res <- results_list[[day]]
    if (myc_gene %in% rownames(res)) {
      myc_sig <- rbind(myc_sig, data.frame(
        Age = age,
        Day = as.numeric(gsub("Day", "", day)),
        FDR = res[myc_gene, "adj.P.Val"],
        significant = res[myc_gene, "adj.P.Val"] < 0.05
      ))
    }
  }
}

myc_summary <- myc_summary %>%
  left_join(myc_sig, by = c("Age", "Day"))

fig5a <- ggplot(myc_summary, aes(x = Day, y = mean_logFC, color = Age, fill = Age, group = Age)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = mean_logFC - se_logFC, ymax = mean_logFC + se_logFC), 
              alpha = 0.3, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(data = myc_summary %>% filter(significant), size = 3.5, stroke = 1) +
  geom_point(data = myc_summary %>% filter(!significant),
             shape = 21, size = 3.5, stroke = 1, fill = "white") +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  scale_x_continuous(breaks = c(2, 10, 20, 30)) +
  labs(
    x = "Day",
    y = "log2FC (Stim vs Control)",
    title = "A: MYC Expression",
    subtitle = "Ribbon = SE, filled = FDR < 0.05"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ============================================================
# PANEL B: MYC TARGETS V1 AND V2
# ============================================================

gsea_myc <- gsea_all %>%
  filter(pathway %in% c("MYC TARGETS V1", "MYC TARGETS V2")) %>%
  mutate(
    Day_num = as.numeric(gsub("Day", "", Day)),
    significant = padj < 0.05
  )

fig5b <- ggplot(gsea_myc, aes(x = Day_num, y = NES, color = Age, fill = Age, group = Age)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  geom_point(data = gsea_myc %>% filter(significant), shape = 21, size = 3.5, stroke = 1) +
  geom_point(data = gsea_myc %>% filter(!significant),
             shape = 21, size = 3.5, stroke = 1, fill = "white") +
  facet_wrap(~ pathway, ncol = 2) +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  scale_x_continuous(breaks = c(2, 10, 20, 30)) +
  labs(
    x = "Day",
    y = "NES",
    title = "B: MYC Target Pathway Enrichment",
    subtitle = "Filled = FDR < 0.05, Open = FDR >= 0.05"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

    # ============================================================
    # PANEL C: TOP DIVERGENT MYC TARGETS
    # ============================================================
    
    # Functional category assignment (4 categories)
    assign_functional_category <- function(gene) {
      gene_upper <- toupper(gene)
      
      # Ribosome biogenesis (includes RNA processing, translation machinery, chaperones)
      if (grepl("^RPL[0-9]|^RPS[0-9]|^MRPL|^MRPS|^NOP[0-9]|^RRP|^DDX|^DHX|^LSM[0-9]|^SNRP|^UTP|^IMP[0-9]|^NHP2|^GAR1|^FBL$|^DKC1|^NOP56|^NOP58|^PPAN|^NOC|^EMG1|^WDR|^BOP1|^PES1|^EIF|^EEF|^PAIP|^RPP|^CCT[0-9]|^TCP1|^HSP[ABDE]|^HSPA|^HSPB|^HSPD|^HSPE|^HSPH|^DNAJ|^PPID|^PPIA$|^PPIB$|^PPIF|^FKBP|^CALR|^CANX|^PDIA|^ERP[0-9]", gene_upper)) {
        return("Ribosome biogenesis")
      }
      
      # Cell proliferation (cell cycle, DNA replication, nucleotide synthesis)
      if (grepl("^MCM[0-9]|^CDC[0-9]|^CDK[0-9]|^CCNA|^CCNB|^CCND|^CCNE|^E2F[0-9]|^ORC[0-9]|^MAD2|^BUB|^PCNA|^CDT1|^DBF4|^PLK|^AURK|^TOP2|^RFC[0-9]|^POLA|^POLD|^POLE|^PRIM|^GINS|^FEN1|^LIG1|^RRM[0-9]|^PRPS[0-9]|^IMPDH|^CTPS|^GMPS|^TYMS|^ODC1|^DCTPP|^CAD$|^UMPS|^ATIC|^GART|^PAICS|^PFAS|^ADSL|^ADSS|^PPAT|^NME[0-9]|^TK1$|^DTYMK|^UCK|^DCK$|^UNG", gene_upper)) {
        return("Cell proliferation")
      }
      
      # Metabolism
      if (grepl("^SLC[0-9]|^HK[0-9]|^LDHA|^LDHB|^ENO[0-9]|^PKM|^PGK[0-9]|^ALDOA|^ALDOB|^ALDOC|^TPI1|^GAPDH|^PFK|^GPI$|^PFKP|^PFKM|^PFKL|^PDH|^IDH|^OGDH|^SDHA|^SDHB|^MDH|^CS$|^ACO[0-9]|^FH$", gene_upper)) {
        return("Metabolism")
      }
      
      # Other (everything else)
      return("Other")
    }
    
    # Define category colours (4 categories)
    category_colors <- c(
      "Ribosome biogenesis" = "#2ca02c",
      "Cell proliferation" = "#ff7f0e",
      "Metabolism" = "#ffdd57",
      "Other" = "#d3d3d3"
    )
    
    # Categorise all MYC target genes
    myc_gene_category <- data.frame(
      gene = unique(c(myc_v1_genes, myc_v2_genes))
    ) %>%
      mutate(
        in_v1 = gene %in% myc_v1_genes,
        in_v2 = gene %in% myc_v2_genes,
        functional_category = sapply(gene, assign_functional_category),
        functional_category = factor(functional_category, 
                                     levels = c("Ribosome biogenesis",
                                                "Cell proliferation",
                                                "Metabolism",
                                                "Other"))
      )
    
    cat("\nFunctional categories:\n")
    print(table(myc_gene_category$functional_category))
    
    # Get MYC target genes present in data
    all_myc_genes <- unique(c(myc_v1_genes, myc_v2_genes))
    myc_genes_in_data <- intersect(all_myc_genes, rownames(log2FC_mat))
    
    cat("\nMYC target genes in data:", length(myc_genes_in_data), "\n")
    
    # Calculate divergence for each gene
    myc_divergence <- data.frame()
    
    for (gene in myc_genes_in_data) {
      gene_data <- data.frame()
      
      for (age in c("Young", "Old")) {
        results_list <- if (age == "Young") results_young else results_old
        for (day in c("Day2", "Day10", "Day20", "Day30")) {
          res <- results_list[[day]]
          if (gene %in% rownames(res)) {
            gene_data <- rbind(gene_data, data.frame(
              gene = gene,
              Age = age,
              Day = day,
              Day_num = as.numeric(gsub("Day", "", day)),
              logFC = res[gene, "logFC"],
              FDR = res[gene, "adj.P.Val"]
            ))
          }
        }
      }
      
      if (nrow(gene_data) == 8) {
        gene_wide <- gene_data %>%
          dplyr::select(gene, Day, Age, logFC) %>%
          pivot_wider(names_from = Age, values_from = logFC) %>%
          mutate(diff = Young - Old)
        
        # Max divergence across all timepoints
        max_div <- max(abs(gene_wide$diff), na.rm = TRUE)
        
        myc_divergence <- rbind(myc_divergence, data.frame(
          gene = gene,
          max_divergence = max_div
        ))
      }
    }
    
    # Get top 36 divergent genes (by max divergence)
    top_myc_genes <- myc_divergence %>%
      arrange(desc(max_divergence)) %>%
      head(36) %>%
      pull(gene)
    
    cat("Top 36 divergent MYC targets selected\n")
    
    # Get individual animal data for plotting
    myc_individual_genes <- data.frame()
    
    for (gene in top_myc_genes) {
      if (gene %in% rownames(log2FC_mat)) {
        for (aid in colnames(log2FC_mat)) {
          animal_meta <- mint_meta %>% filter(animal_id == aid)
          if (nrow(animal_meta) == 1) {
            myc_individual_genes <- rbind(myc_individual_genes, data.frame(
              gene = gene,
              animal_id = aid,
              Age = as.character(animal_meta$Age_group),
              Day_num = as.numeric(as.character(animal_meta$time)),
              logFC = log2FC_mat[gene, aid]
            ))
          }
        }
      }
    }
    
    # Add functional category
    myc_individual_genes <- myc_individual_genes %>%
      left_join(myc_gene_category %>% dplyr::select(gene, functional_category), by = "gene")
    
    # Get significance info
    sig_info <- data.frame()
    for (gene in top_myc_genes) {
      for (age in c("Young", "Old")) {
        results_list <- if (age == "Young") results_young else results_old
        for (day in c("Day2", "Day10", "Day20", "Day30")) {
          res <- results_list[[day]]
          if (gene %in% rownames(res)) {
            sig_info <- rbind(sig_info, data.frame(
              gene = gene,
              Age = age,
              Day_num = as.numeric(gsub("Day", "", day)),
              significant = res[gene, "adj.P.Val"] < 0.05
            ))
          }
        }
      }
    }
    
    # Order genes by max divergence (most to least divergent)
    gene_cat_ordered <- myc_gene_category %>%
      filter(gene %in% top_myc_genes) %>%
      left_join(myc_divergence, by = "gene") %>%
      arrange(desc(max_divergence))
    
    gene_order_by_div <- gene_cat_ordered$gene
    
    cat("\nTop 36 divergent genes by functional category:\n")
    print(table(gene_cat_ordered$functional_category))
    
    # Calculate summary with gene order by divergence
    myc_gene_summary <- myc_individual_genes %>%
      group_by(gene, Age, Day_num, functional_category) %>%
      summarise(
        mean_logFC = mean(logFC),
        se_logFC = sd(logFC) / sqrt(n()),
        n = n(),
        .groups = "drop"
      ) %>%
      left_join(sig_info, by = c("gene", "Age", "Day_num")) %>%
      mutate(gene = factor(gene, levels = gene_order_by_div))
    
    # Create strip background colours - one per gene in order
    strip_fills <- category_colors[as.character(gene_cat_ordered$functional_category)]
    
    # Build the plot first without custom strips
    fig5c_base <- ggplot(myc_gene_summary, aes(x = Day_num, y = mean_logFC, color = Age, fill = Age, group = Age)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_ribbon(aes(ymin = mean_logFC - se_logFC, ymax = mean_logFC + se_logFC),
                  alpha = 0.3, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_point(data = myc_gene_summary %>% filter(significant),
                 shape = 21, size = 2, stroke = 0.8) +
      geom_point(data = myc_gene_summary %>% filter(!significant),
                 shape = 21, size = 2, stroke = 0.8, fill = "white") +
      facet_wrap(~ gene, ncol = 6, scales = "free_y") +
      scale_color_manual(values = age_colors) +
      scale_fill_manual(values = age_colors) +
      scale_x_continuous(breaks = c(2, 10, 20, 30)) +
      labs(
        x = "Day",
        y = "log2FC"
      ) +
      theme_classic(base_size = 9) +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "gray90", colour = "gray50"),
        strip.text = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 7),
        panel.spacing = unit(0.3, "lines")
      )
    
    # Convert to gtable and modify strip colours
    g <- ggplot_gtable(ggplot_build(fig5c_base))
    
    # Find strip-t elements (top strips from facet_wrap)
    strip_both <- which(grepl("strip-t", g$layout$name))
    
    # Apply colours to each strip
    for (i in seq_along(strip_both)) {
      k <- strip_both[i]
      g$grobs[[k]]$grobs[[1]]$children[[1]]$gp$fill <- strip_fills[i]
    }
    
    # Create category legend as separate plot
    legend_data <- data.frame(
      category = factor(names(category_colors), levels = names(category_colors)),
      x = 1:4,
      y = 1
    )
    
    category_legend <- ggplot(legend_data, aes(x = x, y = y, fill = category)) +
      geom_tile(width = 0.9, height = 0.4, color = "white") +
      scale_fill_manual(values = category_colors, name = "") +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8)
      ) +
      guides(fill = guide_legend(nrow = 1))
    
    # Extract just the legend
    legend_grob <- cowplot::get_legend(category_legend)
    
    # Create title
    title_grob <- textGrob(
      "C: Top Divergent MYC Targets (ordered by max divergence)",
      gp = gpar(fontsize = 11, fontface = "bold"),
      hjust = 0, x = 0.02
    )
    
    # Combine: title + plot + legend
    fig5c <- arrangeGrob(
      title_grob,
      g,
      legend_grob,
      ncol = 1,
      heights = c(0.5, 10, 0.8)
    )

# ============================================================
# COMBINE AND SAVE FIGURE 5
# ============================================================

fig5_top <- (fig5a | fig5b) + plot_layout(widths = c(1, 1.5))

# Wrap the gtable for patchwork
fig5c_wrapped <- wrap_elements(full = fig5c)

fig5 <- fig5_top / fig5c_wrapped + plot_layout(heights = c(1, 2.2))

ggsave("paper_v4/figures/Figure5_MYC_execution.png", fig5, width = 14, height = 18, dpi = 300, bg = "white")

cat("\n=== FIGURE 5 COMPLETE ===\n")

# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(myc_gene_category, "paper_v4/results/myc_gene_categories.csv", row.names = FALSE)
write.csv(myc_divergence, "paper_v4/results/myc_gene_divergence.csv", row.names = FALSE)
write.csv(gsea_myc %>% dplyr::select(-leadingEdge), "paper_v4/results/gsea_myc_v1_v2.csv", row.names = FALSE)

# Print category breakdown
cat("\nTop 36 divergent genes by functional category:\n")
print(table(gene_cat_ordered$functional_category))

# ============================================================
# FIGURE 6: SHAVLAKADZE VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 6: SHAVLAKADZE VALIDATION\n")
cat("============================================================\n")

# ============================================================
# PREPARE SHAVLAKADZE DATA
# ============================================================

# Filter to only 6Mo and 24Mo for direct comparison
# Set Young as reference so coefficient = Old - Young
shav_meta_filtered <- meta_shav %>%
  filter(age_months %in% c("6Mo", "24Mo")) %>%
  mutate(Age_group = factor(ifelse(age_months == "6Mo", "Young", "Old"), 
                            levels = c("Young", "Old")))

# Normalise Shavlakadze counts
shav_samples <- shav_meta_filtered$sample_id
shav_counts_subset <- counts_shav[, shav_samples]

dge_shav <- calcNormFactors(DGEList(counts = shav_counts_subset))
shav_logcpm <- cpm(dge_shav, log = TRUE, prior.count = 2)

cat("Shavlakadze samples by age (6Mo vs 24Mo):\n")
print(table(shav_meta_filtered$Age_group))

# ============================================================
# DIFFERENTIAL EXPRESSION: OLD VS YOUNG AT BASELINE
# ============================================================

design_shav <- model.matrix(~ Age_group, data = shav_meta_filtered)
cat("\nDesign matrix columns:", colnames(design_shav), "\n")

fit_shav <- lmFit(shav_logcpm, design_shav)
fit_shav <- eBayes(fit_shav)
shav_de_results <- topTable(fit_shav, coef = 2, number = Inf, sort.by = "none")

n_up <- sum(shav_de_results$adj.P.Val < 0.05 & shav_de_results$logFC > 0)
n_down <- sum(shav_de_results$adj.P.Val < 0.05 & shav_de_results$logFC < 0)
n_total_sig <- n_up + n_down

cat("\nShavlakadze DE results (Old vs Young baseline):\n")
cat("  Up in Old:", n_up, "\n")
cat("  Down in Old:", n_down, "\n")
cat("  Total DEGs:", n_total_sig, "\n")

# ============================================================
# PANEL A: VOLCANO PLOT
# ============================================================

volcano_df <- shav_de_results %>%
  rownames_to_column("gene") %>%
  mutate(
    sig = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up in Old",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down in Old",
      TRUE ~ "NS"
    ),
    neg_log10p = -log10(adj.P.Val)
  )

fig6a <- ggplot(volcano_df, aes(x = logFC, y = neg_log10p, color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Up in Old" = "#D73027", "Down in Old" = "#4575B4", "NS" = "gray80")) +
  scale_x_continuous(limits = c(-5, 5)) +
  annotate("text", x = 3, y = max(volcano_df$neg_log10p) * 0.95, 
           label = paste0("Up in Old: ", n_up), color = "#D73027", size = 4, hjust = 0) +
  annotate("text", x = -3, y = max(volcano_df$neg_log10p) * 0.95, 
           label = paste0("Down in Old: ", n_down), color = "#4575B4", size = 4, hjust = 1) +
  labs(
    x = "log2FC (Old / Young)", 
    y = "-log10(FDR)", 
    title = "A: Baseline Age Differences (Shavlakadze 2019)",
    subtitle = paste0("6Mo vs 24Mo sedentary muscle, ", n_total_sig, " DEGs (FDR < 0.05)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none", 
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ============================================================
# PANEL B: GSEA - SAME PATHWAYS AS FIGURE 3
# ============================================================

shav_ranked <- shav_de_results$t
names(shav_ranked) <- rownames(shav_de_results)
shav_ranked <- sort(shav_ranked[!is.na(shav_ranked)], decreasing = TRUE)

set.seed(42)
shav_gsea <- as.data.frame(fgsea(hallmark_list, shav_ranked, minSize = 15, maxSize = 500, nPermSimple = 10000))

shav_gsea <- shav_gsea %>%
  mutate(
    pathway = gsub("HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway)
  ) %>%
  arrange(pval)

# Filter to Figure 3 pathways
shav_gsea_fig3 <- shav_gsea %>%
  filter(pathway %in% top15_divergent) %>%
  mutate(
    pathway = factor(pathway, levels = rev(top15_divergent)),
    Direction = ifelse(NES > 0, "Up in Old", "Down in Old"),
    significant = padj < 0.05
  )

cat("\nShavlakadze GSEA for Figure 3 pathways:\n")
print(shav_gsea_fig3[, c("pathway", "NES", "padj")])

fig6b <- ggplot(shav_gsea_fig3, aes(x = NES, y = pathway, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label = ifelse(significant, "*", "")), 
            x = ifelse(shav_gsea_fig3$NES > 0, shav_gsea_fig3$NES + 0.1, shav_gsea_fig3$NES - 0.1),
            size = 5, vjust = 0.75) +
  scale_fill_manual(values = c("Up in Old" = "#D73027", "Down in Old" = "#4575B4")) +
  labs(
    x = "NES", 
    y = "", 
    title = "B: Baseline Enrichment (Figure 3 Pathways)",
    subtitle = "* = FDR < 0.05"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none", 
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ============================================================
# PANEL C: BASELINE PREDICTS RESPONSE DIVERGENCE (FIGURE 3 PATHWAYS ONLY)
# ============================================================

# Get mean response divergence for Figure 3 pathways only
prt4_response_avg <- hallmark_combined %>%
  filter(pathway %in% top15_divergent) %>%
  group_by(pathway) %>%
  summarise(
    NES_diff_response = mean(NES_old - NES_young, na.rm = TRUE),
    NES_diff_se = sd(NES_old - NES_young, na.rm = TRUE) / sqrt(n()),
    n_sig_timepoints = sum(padj_young < 0.05 | padj_old < 0.05),
    .groups = "drop"
  )

# Shavlakadze baseline
shav_baseline <- shav_gsea %>%
  dplyr::select(pathway, NES, padj) %>%
  rename(NES_baseline = NES, padj_baseline = padj)

# Combine - only Figure 3 pathways
comparison_df <- prt4_response_avg %>%
  inner_join(shav_baseline, by = "pathway")

cat("\nFigure 3 pathways with data in both datasets:", nrow(comparison_df), "\n")

# Shorten pathway names
comparison_df <- comparison_df %>%
  mutate(
    pathway_short = case_when(
      pathway == "OXIDATIVE PHOSPHORYLATION" ~ "OXPHOS",
      pathway == "MYC TARGETS V1" ~ "MYC V1",
      pathway == "MYC TARGETS V2" ~ "MYC V2",
      pathway == "FATTY ACID METABOLISM" ~ "FA Met",
      pathway == "PI3K AKT MTOR SIGNALING" ~ "PI3K/mTOR",
      pathway == "WNT BETA CATENIN SIGNALING" ~ "Wnt",
      pathway == "XENOBIOTIC METABOLISM" ~ "Xenobiotic",
      pathway == "HEME METABOLISM" ~ "Heme",
      pathway == "G2M CHECKPOINT" ~ "G2M",
      pathway == "E2F TARGETS" ~ "E2F",
      pathway == "KRAS SIGNALING DN" ~ "KRAS DN",
      pathway == "PROTEIN SECRETION" ~ "Secretion",
      pathway == "DNA REPAIR" ~ "DNA Repair",
      pathway == "NOTCH SIGNALING" ~ "Notch",
      pathway == "ADIPOGENESIS" ~ "Adipogenesis",
      pathway == "PEROXISOME" ~ "Peroxisome",
      TRUE ~ pathway
    )
  )

# Calculate correlation with only Figure 3 pathways
cor_test <- cor.test(comparison_df$NES_baseline, comparison_df$NES_diff_response)

cat("\nCorrelation (Figure 3 pathways only):\n")
cat("  n pathways =", nrow(comparison_df), "\n")
cat("  r =", round(cor_test$estimate, 3), "\n")
cat("  p =", signif(cor_test$p.value, 3), "\n")

# Print the data
cat("\nPathway data:\n")
print(comparison_df %>% dplyr::select(pathway_short, NES_baseline, NES_diff_response) %>% arrange(NES_baseline))

fig6c <- ggplot(comparison_df, aes(x = NES_baseline, y = NES_diff_response)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_smooth(method = "lm", se = TRUE, color = "#D73027", fill = "#D73027", linetype = "dashed", alpha = 0.2) +
  geom_point(shape = 21, size = 4, fill = "black", color = "black", alpha = 0.8) +
  geom_text_repel(
    aes(label = pathway_short), 
    size = 3, 
    max.overlaps = 20,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "gray50",
    segment.size = 0.3,
    min.segment.length = 0,
    force = 2,
    force_pull = 0.5,
    direction = "both",
    seed = 42
  ) +
  annotate("text", 
           x = 2.5, 
           y = 1.8,  # Changed from -1.8 to 1.8
           label = paste0("r = ", round(cor_test$estimate, 2), "\np = ", signif(cor_test$p.value, 2)),
           hjust = 1, vjust = 1,  # Changed vjust from 0 to 1
           size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(-3.5, 3.5)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    x = "Baseline NES (Old/Young)", 
    y = "Response dNES (Old - Young)",
    title = "C: Baseline State Predicts Response Divergence",
    subtitle = "Shavlakadze baseline vs PRT4 response (n = 15 pathways)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "none"
  )

# ============================================================
# PANEL D: RESPONSE SIGNATURE AT BASELINE
# ============================================================

comp1_genes <- comp1_loadings$gene
comp1_weights <- setNames(comp1_loadings$loading, comp1_loadings$gene)
genes_in_shav <- intersect(comp1_genes, rownames(shav_logcpm))

cat("\nsPLS-DA Component 1 validation:\n")
cat("  Total Comp1 genes:", length(comp1_genes), "\n")
cat("  Genes in Shavlakadze:", length(genes_in_shav), "\n")

if (length(genes_in_shav) > 50) {
  
  # Calculate weighted score
  expr_mat <- shav_logcpm[genes_in_shav, ]
  expr_z <- t(scale(t(expr_mat)))
  loadings_subset <- comp1_weights[genes_in_shav]
  comp1_scores <- as.numeric(t(loadings_subset) %*% expr_z)
  names(comp1_scores) <- colnames(expr_z)
  
  score_df <- data.frame(
    sample_id = names(comp1_scores),
    Comp1_score = comp1_scores
  ) %>%
    left_join(shav_meta_filtered[, c("sample_id", "Age_group")], by = "sample_id")
  
  comp1_test <- t.test(Comp1_score ~ Age_group, data = score_df)
  
  cat("  Component 1 score t-test: p =", signif(comp1_test$p.value, 3), "\n")
  cat("  Young mean:", round(mean(score_df$Comp1_score[score_df$Age_group == "Young"]), 2), "\n")
  cat("  Old mean:", round(mean(score_df$Comp1_score[score_df$Age_group == "Old"]), 2), "\n")
  
  fig6d <- ggplot(score_df, aes(x = Age_group, y = Comp1_score, fill = Age_group)) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Young" = "#4575B4", "Old" = "#D73027")) +
    annotate("text", x = 1.5, y = max(score_df$Comp1_score) + 0.5,
             label = paste0("p = ", signif(comp1_test$p.value, 2)),
             size = 4, fontface = "bold") +
    labs(
      x = "", 
      y = "Response Signature Score",
      title = "D: Baseline Expression of Response Genes",
      subtitle = "Positive = Young-response genes; Negative = Old-response genes"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "gray40")
    )
  
  write.csv(score_df, "paper_v4/results/shavlakadze_comp1_scores.csv", row.names = FALSE)
  
} else {
  fig6d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient gene overlap", size = 5) +
    theme_void() +
    ggtitle("D: Response Signature at Baseline")
}

# ============================================================
# COMBINE AND SAVE FIGURE 6
# ============================================================

fig6_top <- (fig6a | fig6b) + plot_layout(widths = c(1, 1))
fig6_bottom <- (fig6c | fig6d) + plot_layout(widths = c(1.2, 1))
fig6 <- fig6_top / fig6_bottom + plot_layout(heights = c(1, 1))

ggsave("paper_v4/figures/Figure6_shavlakadze_validation.png", fig6, width = 14, height = 12, dpi = 300, bg = "white")

cat("\n=== FIGURE 6 COMPLETE ===\n")

# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(comparison_df, "paper_v4/results/shavlakadze_baseline_response_comparison.csv", row.names = FALSE)
write.csv(shav_gsea %>% dplyr::select(-leadingEdge), "paper_v4/results/shavlakadze_gsea_results.csv", row.names = FALSE)
write.csv(volcano_df, "paper_v4/results/shavlakadze_de_results.csv", row.names = FALSE)

