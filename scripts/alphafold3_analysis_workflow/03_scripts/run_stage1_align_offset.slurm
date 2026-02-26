#!/bin/bash
#SBATCH --job-name=af3_stage1_align
#SBATCH --output=../04_logs/af3_align_%A_%a.out
#SBATCH --error=../04_logs/af3_align_%A_%a.err
#
# --- The array size will be provided on the command line ---
#
# --- Partition settings from cluster example ---
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00         # <-- INCREASED to 12h to prevent time-limit cancellations
#SBATCH --mem=32gb

# --- OFFSET LOGIC ---
# If LINE_OFFSET is not set by the sbatch command, default it to 0 (for the first batch)
: "${LINE_OFFSET:=0}"

# Calculate the true line number to process
# SLURM_ARRAY_TASK_ID will be 1, 2, 3...
# LINE_OFFSET will be 0, 1000, 2000...
# Example: Task 1, Offset 1000 -> LINE_NUM = 1 + 1000 = 1001
LINE_NUM=$((SLURM_ARRAY_TASK_ID + LINE_OFFSET))

echo "Starting AF3 Stage-1 (Alignment) task $SLURM_ARRAY_TASK_ID"
echo "LINE_OFFSET is $LINE_OFFSET"
echo "Processing line $LINE_NUM of input_list.txt on $(hostname)"
echo "Job started at: $(date)"

# --- 1. Set up paths ---
BASE_DIR=$(pwd)/..
INPUT_LIST_FILE="$BASE_DIR/03_scripts/input_list.txt"
JSON_INPUT_DIR="$BASE_DIR/02_inputs_json"

# Set environment variable for your personal model parameters directory
ALPHAFOLD_MODEL_DIR=/home/hd/hd_hd/hd_bc242/af3-models

# --- ROBUSTNESS FIX: Hard-code the workspace path ---
# ALPHAFOLD_RESULTS_DIR="/gpfs/bwfor/work/ws/hd_bc242-alphafold"
ALPHAFOLD_RESULTS_DIR="/home/hd/hd_hd/hd_bc242/hd_bc242-alphafold"

if [ -z "$ALPHAFOLD_RESULTS_DIR" ]; then
    echo "Error: Hard-coded workspace path is empty?"
    exit 1
fi

# --- 2. Get the specific input file for this task ---
# Use our new LINE_NUM variable
INPUT_BASENAME=$(sed -n "${LINE_NUM}p" $INPUT_LIST_FILE)
if [ -z "$INPUT_BASENAME" ]; then
    echo "Error: No input basename found for task ID $SLURM_ARRAY_TASK_ID (Line $LINE_NUM)"
    exit 1
fi

# Define paths for this task
JSON_INPUT_PATH="$JSON_INPUT_DIR/${INPUT_BASENAME}.json"
OUTPUT_DIR=$ALPHAFOLD_RESULTS_DIR

echo "Input JSON: $JSON_INPUT_PATH"
echo "Output base directory: $OUTPUT_DIR"
echo "Model parameters directory: $ALPHAFOLD_MODEL_DIR"

# --- 3. Load software module ---
module load bio/alphafold/3.0.1

# --- 4. Run Stage 1 (Alignment) ---
python $ALPHAFOLD_BIN_DIR/run_alphafold.py \
    --json_path=$JSON_INPUT_PATH \
    --db_dir=$ALPHAFOLD_DATABASES \
    --model_dir=$ALPHAFOLD_MODEL_DIR \
    --output_dir=$OUTPUT_DIR \
    --norun_inference

echo "Stage 1 finished at: $(date)"
