#!/bin/bash
#SBATCH -J hdmanova
#SBATCH -p barbun
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module purge
module load apps/R/4.5.2-gcc-14.1.0

cd "${SLURM_SUBMIT_DIR}"

GRID_FILE="${GRID_FILE:?GRID_FILE not set}"
STUDY="${STUDY:?STUDY not set}"
OUTDIR="${OUTDIR:?OUTDIR not set}"

BASE_SEED="${BASE_SEED:-12345}"
R_REPS="${R_REPS:-1000}"
R_BLOCK="${R_BLOCK:-100}"
B_PERM="${B_PERM:-199}"
ALPHA_LEVEL="${ALPHA_LEVEL:-0.05}"
ETA="${ETA:-0.75}"
DELTA_ROB="${DELTA_ROB:-4.0}"
RHO="${RHO:-0.5}"
MC_CORES="${MC_CORES:-20}"

TASK_ID="${HD_TASK_ID:-$SLURM_ARRAY_TASK_ID}"

echo "PWD=$(pwd)"
echo "GRID_FILE=$GRID_FILE"
echo "STUDY=$STUDY"
echo "OUTDIR=$OUTDIR"
echo "TASK_ID=$TASK_ID"
echo "R_REPS=$R_REPS"
echo "R_BLOCK=$R_BLOCK"
echo "MC_CORES=$MC_CORES"

mkdir -p "$OUTDIR"
ls -ld "$OUTDIR"

Rscript code/simulation/run_hdmanova_chunk.R \
  --grid "$GRID_FILE" \
  --study "$STUDY" \
  --task-id "$TASK_ID" \
  --outdir "$OUTDIR" \
  --R "$R_REPS" \
  --R-block "$R_BLOCK" \
  --B "$B_PERM" \
  --alpha "$ALPHA_LEVEL" \
  --seed "$BASE_SEED" \
  --eta "$ETA" \
  --delta-rob "$DELTA_ROB" \
  --rho "$RHO" \
  --mc-cores "$MC_CORES"