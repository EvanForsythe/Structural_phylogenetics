#!/bin/bash

#SBATCH --partition=forsythe.q
#SBATCH --account=forsythe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a40:1
#SBATCH --time=12:00:00
#SBATCH --job-name=colabfold
#SBATCH --mail-user=evan.forsythe@osucascades.edu
#SBATCH --mail-type=END,FAIL


# ============================================================
# Usage
# ============================================================
#
# sbatch 02_colabfold.sbatch <input_fasta> <output_dir>
#
# Example:
#
# sbatch 02_colabfold.sbatch \
#   /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family/HOG0015056.fa \
#   /scratch/forsythe/Structural_phylogenetics/RNAase_T2_RNS1_family
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
    echo "sbatch 02_colabfold.sbatch <input_fasta> <output_dir>"
    exit 1
fi


# ------------------------------------------------------------
# Convert paths to absolute paths
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

COLABFOLD_DIR="${OUTPUT_DIR}/colabfold_results"
LOG_DIR="${OUTPUT_DIR}/logs"

COLABFOLD_DATA="/scratch/forsythe/Structural_phylogenetics/colabfold_cache"

mkdir -p "$COLABFOLD_DIR"
mkdir -p "$LOG_DIR"


# ------------------------------------------------------------
# Redirect stdout/stderr to project-specific log
# ------------------------------------------------------------

LOG_FILE="${LOG_DIR}/colabfold_${SLURM_JOB_ID}.log"

exec > >(tee -a "$LOG_FILE") 2>&1


# ------------------------------------------------------------
# Activate ColabFold environment
# ------------------------------------------------------------

source /opt/ohpc/pub/apps/anaconda/3.12/etc/profile.d/conda.sh

conda activate /scratch/forsythe/conda_envs/colabfold


# ------------------------------------------------------------
# Print job information
# ------------------------------------------------------------

echo
echo "============================================================"
echo "ColabFold structure prediction"
echo "============================================================"
echo "Date:              $(date)"
echo "SLURM job ID:      ${SLURM_JOB_ID}"
echo "Node:              $(hostname)"
echo
echo "Input FASTA:       ${INPUT_FASTA}"
echo "Output directory:  ${OUTPUT_DIR}"
echo "ColabFold output:  ${COLABFOLD_DIR}"
echo "ColabFold data:    ${COLABFOLD_DATA}"
echo "Log file:          ${LOG_FILE}"
echo
echo "ColabFold:"
which colabfold_batch
echo
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
echo

nvidia-smi

echo
echo "============================================================"
echo "Starting ColabFold"
echo "============================================================"
echo


# ------------------------------------------------------------
# Run ColabFold
# ------------------------------------------------------------

colabfold_batch \
    "$INPUT_FASTA" \
    "$COLABFOLD_DIR" \
    --data "$COLABFOLD_DATA"

EXIT_STATUS=$?


# ------------------------------------------------------------
# Report completion
# ------------------------------------------------------------

echo
echo "============================================================"

if [[ $EXIT_STATUS -eq 0 ]]; then
    echo "ColabFold completed successfully."
else
    echo "ERROR: ColabFold exited with status ${EXIT_STATUS}."
fi

echo "Finished: $(date)"
echo "============================================================"
echo

exit $EXIT_STATUS