#!/bin/bash
#SBATCH --job-name=af3_stage2_infer
#SBATCH --output=../04_logs/af3_infer_%A_%a.out
#SBATCH --error=../04_logs/af3_infer_%A_%a.err
#
# --- Partition settings ---
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1,gpumem_per_gpu=80GB
#SBATCH --time=02:00:00
#SBATCH --mem=100gb

# --- OFFSET LOGIC ---
: "${LINE_OFFSET:=0}"
LINE_NUM=$((SLURM_ARRAY_TASK_ID + LINE_OFFSET))

# --- 1. Set up paths ---
BASE_DIR=$(pwd)/..
INPUT_LIST_FILE="$BASE_DIR/03_scripts/input_list.txt"
ALPHAFOLD_MODEL_DIR=$HOME/af3-models
ALPHAFOLD_RESULTS_DIR="/gpfs/bwfor/work/ws/hd_bc242-alphafold"

# --- 2. Get Input Name ---
INPUT_BASENAME=$(sed -n "${LINE_NUM}p" $INPUT_LIST_FILE)
if [ -z "$INPUT_BASENAME" ]; then
    echo "Error: No input basename found for line $LINE_NUM"
    exit 1
fi

echo "Processing: $INPUT_BASENAME (Line $LINE_NUM)"

# --- 3. Robust Path Finding (Flexible Matching) ---

# Convert Input Name to the format seen on disk
# 1. Lowercase everything (e.g., ppi_Hsp90aa1_and_Ttn -> ppi_hsp90aa1_and_ttn)
BASENAME_LOWER=$(echo "$INPUT_BASENAME" | tr '[:upper:]' '[:lower:]')

# 2. Create a version without "_and_" (e.g., ppi_hsp90aa1_ttn)
#    This matches the pattern seen in your ls output (e.g., ppi_aak1_ap2m1_...)
BASENAME_NO_AND=${BASENAME_LOWER//_and_/_}

# 3. Find the directory
# We look for a directory starting with the "no_and" pattern
FOUND_DIR=$(find "$ALPHAFOLD_RESULTS_DIR" -maxdepth 1 -type d -name "${BASENAME_NO_AND}*" | head -n 1)

# Fallback: If not found, try the version WITH "and" just in case
if [ -z "$FOUND_DIR" ]; then
    FOUND_DIR=$(find "$ALPHAFOLD_RESULTS_DIR" -maxdepth 1 -type d -name "${BASENAME_LOWER}*" | head -n 1)
fi

# Check if we found anything
if [ -z "$FOUND_DIR" ]; then
    echo "SKIP: Alignment directory not found for $INPUT_BASENAME"
    echo "      Searched for: ${BASENAME_NO_AND}* OR ${BASENAME_LOWER}*"
    exit 0 # Exit success to skip
fi

echo "Found Alignment Directory: $FOUND_DIR"
TARGET_DIR="$FOUND_DIR"

# --- 4. RESUME LOGIC ---
# Check if output already exists (summary_confidences.json)
if find "$TARGET_DIR" -name "summary_confidences.json" -print -quit | grep -q .; then
    echo "SKIP: Output already exists in $TARGET_DIR (or subdirectory). Job complete."
    exit 0
fi

# --- 5. Find Data File ---
# We don't guess the name anymore. We find the file ending in _data.json inside the target dir.
JSON_INPUT_PATH=$(find "$TARGET_DIR" -maxdepth 1 -name "*_data.json" | head -n 1)

if [ -z "$JSON_INPUT_PATH" ]; then
    echo "SKIP: _data.json file not found inside $TARGET_DIR"
    exit 0
fi

echo "Found Alignment Data: $JSON_INPUT_PATH"

# --- 6. Run AlphaFold 3 ---
echo "Starting Inference..."
module load bio/alphafold/3.0.1

python $ALPHAFOLD_BIN_DIR/run_alphafold.py \
    --json_path=$JSON_INPUT_PATH \
    --db_dir=$ALPHAFOLD_DATABASES \
    --model_dir=$ALPHAFOLD_MODEL_DIR \
    --output_dir=$ALPHAFOLD_RESULTS_DIR \
    --norun_data_pipeline

echo "Stage 2 finished at: $(date)"
