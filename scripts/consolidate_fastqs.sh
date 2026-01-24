#!/bin/bash

# =============================================================================
# FASTQ Consolidation Script for PRT2/4 and Shavlakadze2019
# =============================================================================

set -e  # Exit on any error

BASE_DIR="/Volumes/Samsung_SSD/rna-seq/raw_data/PRT4_paper_extracted_fastqs"

# Create new directory structure
echo "Creating directory structure..."
mkdir -p "${BASE_DIR}/combined_fastqs/PRT2_4"
mkdir -p "${BASE_DIR}/combined_fastqs/Shavlakadze2019"
mkdir -p "${BASE_DIR}/archive"

# =============================================================================
# PART 1: PRT2 - Concatenate lanes and rename
# =============================================================================
echo ""
echo "Processing PRT2 samples (concatenating 4 lanes each)..."

PRT2_DIR="${BASE_DIR}/PRT2_extracted_fastqs"
OUT_DIR="${BASE_DIR}/combined_fastqs/PRT2_4"

# Get unique sample prefixes
for prefix in $(ls "${PRT2_DIR}" | grep "_R1_" | cut -d'_' -f1 | sort -u); do
    
    # Extract animal number and side from prefix
    # Handle the 265-L typo case
    if [[ "$prefix" == "GC-JJ-9173-265-L-17" ]]; then
        animal="265"
        side="L"
    else
        # Standard format: GC-JJ-9173-257L-1
        # Extract the number+side portion (e.g., 257L)
        numsid=$(echo "$prefix" | sed 's/GC-JJ-9173-//' | sed 's/-[0-9]*$//')
        animal=$(echo "$numsid" | sed 's/[LR]$//')
        side=$(echo "$numsid" | grep -o '[LR]$')
    fi
    
    # Construct new name
    newname="TA${side}${animal}"
    
    echo "  ${prefix} -> ${newname}"
    
    # Concatenate R1 lanes (L001-L004)
    cat "${PRT2_DIR}/${prefix}"_*_R1_*.fastq.gz > "${OUT_DIR}/${newname}_R1.fastq.gz"
    
    # Concatenate R2 lanes (L001-L004)
    cat "${PRT2_DIR}/${prefix}"_*_R2_*.fastq.gz > "${OUT_DIR}/${newname}_R2.fastq.gz"
    
done

echo "PRT2 complete: $(ls ${OUT_DIR} | wc -l) files created"

# =============================================================================
# PART 2: PRT4 Part1 - Rename and move
# =============================================================================
echo ""
echo "Processing PRT4 Part1 samples..."

PRT4_P1_DIR="${BASE_DIR}/PRT4_part1_extracted_fastqs/01.RawData"

for sample_dir in "${PRT4_P1_DIR}"/TAR*; do
    sample=$(basename "$sample_dir")
    
    # Find the fastq files (format: TAR448_EKRO..._L4_1.fq.gz)
    r1=$(ls "${sample_dir}"/*_1.fq.gz 2>/dev/null | grep -v "^\._" | head -1)
    r2=$(ls "${sample_dir}"/*_2.fq.gz 2>/dev/null | grep -v "^\._" | head -1)
    
    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "  ${sample} -> ${sample}_R1/R2.fastq.gz"
        mv "$r1" "${OUT_DIR}/${sample}_R1.fastq.gz"
        mv "$r2" "${OUT_DIR}/${sample}_R2.fastq.gz"
    else
        echo "  WARNING: Missing files for ${sample}"
    fi
done

echo "PRT4 Part1 complete"

# =============================================================================
# PART 3: PRT4 Part2 - Rename and move
# =============================================================================
echo ""
echo "Processing PRT4 Part2 samples..."

PRT4_P2_DIR="${BASE_DIR}/PRT4_part2_extracted_fastqs/01.RawData"

for sample_dir in "${PRT4_P2_DIR}"/*/; do
    sample=$(basename "$sample_dir")
    
    # Find the fastq files
    r1=$(ls "${sample_dir}"/*_1.fq.gz 2>/dev/null | grep -v "^\._" | head -1)
    r2=$(ls "${sample_dir}"/*_2.fq.gz 2>/dev/null | grep -v "^\._" | head -1)
    
    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "  ${sample} -> ${sample}_R1/R2.fastq.gz"
        mv "$r1" "${OUT_DIR}/${sample}_R1.fastq.gz"
        mv "$r2" "${OUT_DIR}/${sample}_R2.fastq.gz"
    else
        echo "  WARNING: Missing files for ${sample}"
    fi
done

echo "PRT4 Part2 complete"

# =============================================================================
# PART 4: Shavlakadze - Rename _1/_2 to _R1/_R2
# =============================================================================
echo ""
echo "Processing Shavlakadze2019 samples..."

SHAV_IN="${BASE_DIR}/Shavlekadze_fastqs"
SHAV_OUT="${BASE_DIR}/combined_fastqs/Shavlakadze2019"

# Copy metadata files
cp "${SHAV_IN}/SraRunTable.csv" "${SHAV_OUT}/"
cp "${SHAV_IN}/SRR_Acc_List.txt" "${SHAV_OUT}/"

# Rename and move fastq files
for fq in "${SHAV_IN}"/SRR*_1.fastq.gz; do
    if [[ -f "$fq" ]]; then
        base=$(basename "$fq" _1.fastq.gz)
        echo "  ${base}_1/2 -> ${base}_R1/R2"
        mv "$fq" "${SHAV_OUT}/${base}_R1.fastq.gz"
        mv "${SHAV_IN}/${base}_2.fastq.gz" "${SHAV_OUT}/${base}_R2.fastq.gz"
    fi
done

echo "Shavlakadze2019 complete"

# =============================================================================
# PART 5: Move original folders to archive
# =============================================================================
echo ""
echo "Moving original folders to archive..."

mv "${BASE_DIR}/PRT2_extracted_fastqs" "${BASE_DIR}/archive/"
mv "${BASE_DIR}/PRT4_part1_extracted_fastqs" "${BASE_DIR}/archive/"
mv "${BASE_DIR}/PRT4_part2_extracted_fastqs" "${BASE_DIR}/archive/"
mv "${BASE_DIR}/Shavlekadze_fastqs" "${BASE_DIR}/archive/"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================="
echo "COMPLETE"
echo "============================================="
echo "PRT2_4 samples: $(ls ${BASE_DIR}/combined_fastqs/PRT2_4/*.fastq.gz | wc -l) files"
echo "Shavlakadze2019 samples: $(ls ${BASE_DIR}/combined_fastqs/Shavlakadze2019/*.fastq.gz | wc -l) files"
echo ""
echo "Verify with:"
echo "  ls ${BASE_DIR}/combined_fastqs/PRT2_4/ | head -20"
echo "  ls ${BASE_DIR}/combined_fastqs/Shavlakadze2019/ | head -20"
