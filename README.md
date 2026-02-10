# prt4-aged-muscle-rnaseq

Comparative RNA-seq analysis of young versus old rat skeletal muscle in response to 30 days of daily programmed peripheral nerve stimulation. This repository contains the complete analytical pipeline from raw FASTQ processing through to publication-quality figures and biological interpretation.

> **Note:** This is an active research project. Code is provided to demonstrate analytical methods. Data files are not included as the work is unpublished.

## Study Design

- **Species:** Rat (*Rattus norvegicus*)
- **Age groups:** Young (PRT2 cohort) and Old (PRT4 cohort)
- **Intervention:** Daily programmed peripheral nerve stimulation (30 days)
- **Timepoints:** Days 2, 10, 20, 30
- **Design:** Paired within-animal (stimulated vs. control limb)
- **Tissue:** Tibialis anterior (TA) skeletal muscle
- **Sequencing:** Bulk RNA-seq (young and old cohorts sequenced on different platforms)

## Repository Structure

    prt4-aged-muscle-rnaseq/
    ├── scripts/
    │   ├── consolidate_fastqs.sh          # FASTQ lane concatenation across sequencing batches
    │   ├── run_salmon.sh                  # Salmon pseudo-alignment with GC/sequence bias correction
    │   ├── tximport_salmon.R              # Transcript-to-gene aggregation via tximport + biomaRt
    │   └── PRT4_v4_analysis_250210.R      # Complete downstream analysis pipeline (~2,000 lines)
    ├── qc_reports/                        # MultiQC and FastQC quality control outputs
    └── README.md

## Upstream Pipeline

1. **FASTQ Consolidation** (consolidate_fastqs.sh) — Concatenates FASTQ files across sequencing lanes and batches into per-sample files with consistent naming.
2. **Pseudo-alignment** (run_salmon.sh) — Quantifies transcript abundances using Salmon with GC bias and sequence bias correction against the rat transcriptome.
3. **Gene-level Aggregation** (tximport_salmon.R) — Aggregates transcript-level counts to gene-level using tximport, with Ensembl annotation via biomaRt. Merges experimental data with a public validation dataset (Shavlakadze et al., 2019) into a unified count matrix.

## Downstream Analysis Pipeline

The main analysis script (PRT4_v4_analysis_250210.R) contains the complete downstream workflow:

### Data Preparation (Lines 1-238)
- Package loading, metadata setup, and dataset subsetting
- Gene filtering (minimum count thresholds applied independently per age group)
- Integration of external validation dataset (Shavlakadze et al., 2019)
- QC sample swap correction identified through pipeline checks

### Phenotypic Analysis (Lines 240-433)
- Muscle mass and hypertrophy quantification
- Linear mixed models (lme4/lmerTest) with post-hoc comparisons (emmeans)
- Paired stimulated vs. control limb analysis

### Transcriptome Overview (Lines 435-750)
- PCA analysis with consistent axis scaling across age groups
- Differential expression via LIMMA-voom with duplicateCorrelation for paired designs
- Sex included as covariate in the old cohort where required
- DEG counting across timepoints and age groups
- Quadrant correlation plots comparing young vs. old responses
- Venn diagram overlap analysis

### Pathway Divergence (Lines 752-906)
- Gene set enrichment analysis (fgsea) with MSigDB Hallmark gene sets
- Age-by-timepoint GSEA comparison
- Identification of top divergent pathways between young and old

### Multivariate Discrimination (Lines 908-1185)
- Sparse Partial Least Squares Discriminant Analysis (sPLS-DA via MixOmics)
- Log2FC calculation per animal (stimulated vs. control)
- Score plots, loading extraction, and trajectory analysis
- Gene Ontology enrichment (clusterProfiler) of discriminating features
- STRING network integration

### Transcription Factor and Target Gene Analysis (Lines 1187-1654)
- Transcription factor activity inference (DoRothEA/VIPER)
- Functional gene categorisation
- Single-sample pathway scoring (GSVA)

### External Validation (Lines 1656-1973)
- Independent validation against Shavlakadze et al. (2019) public dataset
- Differential expression at baseline (old vs. young)
- Correlation between baseline differences and stimulation response divergence
- Demonstrates that pre-existing transcriptomic differences predict response outcomes

## Key Analytical Decisions

- **Qualitative comparison:** Because young and old cohorts were sequenced on different platforms, age comparisons use direction and significance rather than direct magnitude comparisons.
- **Paired design:** duplicateCorrelation accounts for within-animal pairing (stimulated vs. control limb).
- **External validation:** Findings validated against an independent public dataset to strengthen biological conclusions.

## Dependencies

R packages: dplyr, tidyr, ggplot2, edgeR, limma, pheatmap, mixOmics, fgsea, msigdbr, patchwork, lme4, lmerTest, emmeans, dorothea, viper, GSVA, clusterProfiler, org.Rn.eg.db, STRINGdb, cowplot, ggrepel

Bash tools: Salmon, MultiQC, FastQC

## Author

Jack Edmondson — PhD Candidate, Liverpool John Moores University

Supervisors: Professor Jonathan Jarvis, Professor Jatin Burniston

## Status

Manuscript in preparation. Analysis pipeline complete.
