#!/bin/bash

#SBATCH --partition=forsythe.q
#SBATCH --account=forsythe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=foldmason
#SBATCH --mail-user=evan.forsythe@osucascades.edu
#SBATCH --mail-type=END,FAIL


# ============================================================
# Usage
# ============================================================
#
# sbatch 04_foldmason.sbatch <input_fasta> <output_dir>
#
# Example:
#
# sbatch 04_foldmason.sbatch \
#   /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa \
#   /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
#
# This script expects Script 3 to have created:
#
#   <output_dir>/structures/
#
# containing exactly one PDB for every sequence in the
# original FASTA.
#
# ============================================================


# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------

INPUT_FASTA="$1"
OUTPUT_DIR="$2"

if [[ -z "$INPUT_FASTA" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Missing required arguments."
    echo
    echo "Usage:"
    echo "sbatch 04_foldmason.sbatch <input_fasta> <output_dir>"
    exit 1
fi


# ------------------------------------------------------------
# Convert to absolute paths
# ------------------------------------------------------------

INPUT_FASTA=$(realpath "$INPUT_FASTA")
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")


# ------------------------------------------------------------
# Check input FASTA
# ------------------------------------------------------------

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA does not exist:"
    echo "$INPUT_FASTA"
    exit 1
fi


# ------------------------------------------------------------
# Define directories
# ------------------------------------------------------------

STRUCTURE_DIR="${OUTPUT_DIR}/structures"
FOLDMASON_DIR="${OUTPUT_DIR}/foldmason"
TMP_DIR="${FOLDMASON_DIR}/tmp"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$FOLDMASON_DIR"
mkdir -p "$TMP_DIR"
mkdir -p "$LOG_DIR"


# ------------------------------------------------------------
# Redirect stdout/stderr to project-specific log
# ------------------------------------------------------------

LOG_FILE="${LOG_DIR}/foldmason_${SLURM_JOB_ID}.log"

exec > >(tee -a "$LOG_FILE") 2>&1


# ------------------------------------------------------------
# Activate structural phylogenetics environment
# ------------------------------------------------------------

source /opt/ohpc/pub/apps/anaconda/3.12/etc/profile.d/conda.sh

conda activate structural_phylo


# ------------------------------------------------------------
# Print job information
# ------------------------------------------------------------

echo
echo "============================================================"
echo "FoldMason structural multiple alignment"
echo "============================================================"
echo "Date:              $(date)"
echo "SLURM job ID:      ${SLURM_JOB_ID}"
echo "Node:              $(hostname)"
echo
echo "Input FASTA:       ${INPUT_FASTA}"
echo "Output directory:  ${OUTPUT_DIR}"
echo "Structure dir:     ${STRUCTURE_DIR}"
echo "FoldMason output:  ${FOLDMASON_DIR}"
echo "Temporary dir:     ${TMP_DIR}"
echo "Log file:          ${LOG_FILE}"
echo


# ------------------------------------------------------------
# Check structure directory
# ------------------------------------------------------------

if [[ ! -d "$STRUCTURE_DIR" ]]; then
    echo "ERROR: Structure directory does not exist:"
    echo "$STRUCTURE_DIR"
    echo
    echo "Run Script 3 first."
    exit 1
fi


# ------------------------------------------------------------
# Check that every FASTA ID has exactly one corresponding PDB
# ------------------------------------------------------------

EXPECTED_COUNT=0
ERROR_COUNT=0
PDB_FILES=()

while read -r SEQUENCE_ID; do

    EXPECTED_COUNT=$((EXPECTED_COUNT + 1))

    PDB_FILE="${STRUCTURE_DIR}/${SEQUENCE_ID}.pdb"

    if [[ ! -f "$PDB_FILE" ]]; then
        echo "ERROR: Missing PDB for:"
        echo "  ${SEQUENCE_ID}"
        ERROR_COUNT=$((ERROR_COUNT + 1))
    else
        PDB_FILES+=("$PDB_FILE")
    fi

done < <(grep '^>' "$INPUT_FASTA" | sed 's/^>//' | awk '{print $1}')


# ------------------------------------------------------------
# Check for unexpected/extra PDB files
# ------------------------------------------------------------

shopt -s nullglob
ALL_PDB_FILES=("${STRUCTURE_DIR}"/*.pdb)
shopt -u nullglob

ACTUAL_COUNT=${#ALL_PDB_FILES[@]}
MATCHED_COUNT=${#PDB_FILES[@]}


echo "Expected proteins: ${EXPECTED_COUNT}"
echo "Matched PDBs:      ${MATCHED_COUNT}"
echo "Total PDB files:   ${ACTUAL_COUNT}"
echo


if [[ "$ACTUAL_COUNT" -ne "$EXPECTED_COUNT" ]]; then

    echo "ERROR: Number of PDB files does not match FASTA records."
    echo "Expected: ${EXPECTED_COUNT}"
    echo "Found:    ${ACTUAL_COUNT}"

    ERROR_COUNT=$((ERROR_COUNT + 1))

fi


if [[ "$ERROR_COUNT" -gt 0 ]]; then

    echo
    echo "============================================================"
    echo "VALIDATION FAILED"
    echo "============================================================"
    echo "FoldMason will not be run."
    echo
    echo "Check the structures/ directory and Script 3 output."
    echo "============================================================"

    exit 1

fi


echo "Structure validation passed."
echo


# ------------------------------------------------------------
# Check FoldMason installation
# ------------------------------------------------------------

echo "FoldMason executable:"
which foldmason

echo
echo "FoldMason version:"
foldmason version

echo


# ------------------------------------------------------------
# Determine family name from FASTA filename
# ------------------------------------------------------------

FASTA_BASENAME=$(basename "$INPUT_FASTA")

# Remove common FASTA extensions
FAMILY_NAME="${FASTA_BASENAME%.fasta}"
FAMILY_NAME="${FAMILY_NAME%.fa}"
FAMILY_NAME="${FAMILY_NAME%.faa}"

RESULT_PREFIX="${FOLDMASON_DIR}/${FAMILY_NAME}_foldmason"


echo "Family name:       ${FAMILY_NAME}"
echo "Result prefix:     ${RESULT_PREFIX}"
echo


# ------------------------------------------------------------
# Run FoldMason
# ------------------------------------------------------------

echo "============================================================"
echo "Starting FoldMason"
echo "============================================================"
echo

foldmason easy-msa \
    "${PDB_FILES[@]}" \
    "$RESULT_PREFIX" \
    "$TMP_DIR" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --report-mode 1

EXIT_STATUS=$?


# ------------------------------------------------------------
# Report completion
# ------------------------------------------------------------

echo
echo "============================================================"

if [[ $EXIT_STATUS -eq 0 ]]; then

    echo "FoldMason completed successfully."
    echo
    echo "Output files:"
    echo

    find "$FOLDMASON_DIR" \
        -maxdepth 1 \
        -type f \
        -printf "  %f\n" \
        | sort

else

    echo "ERROR: FoldMason exited with status ${EXIT_STATUS}."

fi

echo
echo "Finished: $(date)"
echo "============================================================"
echo

exit $EXIT_STATUS