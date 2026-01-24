#!/bin/bash

# =============================================================================
# Salmon Quantification Script for PRT2/4 and Shavlakadze2019
# =============================================================================

set -e

BASE_DIR="/Volumes/Samsung_SSD/rna-seq/raw_data/PRT4_paper_extracted_fastqs"
INDEX="/Volumes/Samsung_SSD/rna-seq/pipeline/references/salmon_index_rn7"
THREADS=6  # Leave 2 cores free for system

# Create output directories
mkdir -p "${BASE_DIR}/salmon_output/PRT2_4"
mkdir -p "${BASE_DIR}/salmon_output/Shavlakadze2019"

# =============================================================================
# PRT2_4 Samples
# =============================================================================
echo "============================================="
echo "Processing PRT2_4 samples"
echo "============================================="

FASTQ_DIR="${BASE_DIR}/combined_fastqs/PRT2_4"
OUT_DIR="${BASE_DIR}/salmon_output/PRT2_4"

# Get unique sample names (strip _R1.fastq.gz and _R2.fastq.gz)
for r1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    
    if [[ -f "$r2" ]]; then
        echo ""
        echo "Processing: ${sample}"
        echo "  R1: $(basename $r1)"
        echo "  R2: $(basename $r2)"
        
        salmon quant \
            -i "$INDEX" \
            -l A \
            -1 "$r1" \
            -2 "$r2" \
            -p $THREADS \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o "${OUT_DIR}/${sample}"
        
        echo "  Done: ${sample}"
    else
        echo "WARNING: Missing R2 for ${sample}"
    fi
done

echo ""
echo "PRT2_4 complete: $(ls -d ${OUT_DIR}/*/ 2>/dev/null | wc -l) samples quantified"

# =============================================================================
# Shavlakadze2019 Samples
# =============================================================================
echo ""
echo "============================================="
echo "Processing Shavlakadze2019 samples"
echo "============================================="

FASTQ_DIR="${BASE_DIR}/combined_fastqs/Shavlakadze2019"
OUT_DIR="${BASE_DIR}/salmon_output/Shavlakadze2019"

for r1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    
    if [[ -f "$r2" ]]; then
        echo ""
        echo "Processing: ${sample}"
        
        salmon quant \
            -i "$INDEX" \
            -l A \
            -1 "$r1" \
            -2 "$r2" \
            -p $THREADS \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o "${OUT_DIR}/${sample}"
        
        echo "  Done: ${sample}"
    else
        echo "WARNING: Missing R2 for ${sample}"
    fi
done

echo ""
echo "Shavlakadze2019 complete: $(ls -d ${OUT_DIR}/*/ 2>/dev/null | wc -l) samples quantified"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================="
echo "QUANTIFICATION COMPLETE"
echo "============================================="
echo "PRT2_4: $(ls -d ${BASE_DIR}/salmon_output/PRT2_4/*/ 2>/dev/null | wc -l) samples"
echo "Shavlakadze2019: $(ls -d ${BASE_DIR}/salmon_output/Shavlakadze2019/*/ 2>/dev/null | wc -l) samples"
echo ""
echo "Output location: ${BASE_DIR}/salmon_output/"
