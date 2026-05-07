#!/bin/bash
set -euo pipefail

# ============================================================
# submit_robustness.sh
# Submit only the Section 3.3 robustness simulation.
# Uses the existing scripts/slurm_hdmanova_array.sbatch.
# No new SLURM file is required.
#
# Usage:
#   bash scripts/submit_robustness.sh pilot
#   bash scripts/submit_robustness.sh main
#
# Optional overrides:
#   DELTA_ROB=8 bash scripts/submit_robustness.sh pilot
#   R_REPS=500 B_PERM=199 bash scripts/submit_robustness.sh pilot
# ============================================================

MODE="${1:-pilot}"

module purge
module load apps/R/4.5.2-gcc-14.1.0

cd "$(dirname "$0")/.."

GRID="results/tables/grid_robustness.csv"
if [ ! -f "$GRID" ]; then
  echo "ERROR: $GRID not found. Put the existing grid_robustness.csv under results/tables." >&2
  exit 1
fi

mkdir -p logs
mkdir -p results/raw_chunks/robustness
mkdir -p results/tables
mkdir -p results/figures/robustness

if [ "$MODE" = "pilot" ]; then
  R_REPS="${R_REPS:-200}"
  R_BLOCK="${R_BLOCK:-100}"
  B_PERM="${B_PERM:-99}"
  DELTA_ROB="${DELTA_ROB:-6.0}"
  OUTDIR="${OUTDIR:-results/raw_chunks/robustness_pilot_delta${DELTA_ROB//./}}"
elif [ "$MODE" = "main" ]; then
  R_REPS="${R_REPS:-1000}"
  R_BLOCK="${R_BLOCK:-100}"
  B_PERM="${B_PERM:-199}"
  DELTA_ROB="${DELTA_ROB:-6.0}"
  OUTDIR="${OUTDIR:-results/raw_chunks/robustness}"
else
  echo "ERROR: MODE must be 'pilot' or 'main'." >&2
  exit 1
fi

NCELLS=$(($(wc -l < "$GRID") - 1))
MAX_ARRAY=99
BASE_SEED="${BASE_SEED:-20260503}"
ALPHA_LEVEL="${ALPHA_LEVEL:-0.05}"
ETA="${ETA:-0.75}"       # ignored for robustness; run_hdmanova_chunk.R uses adaptive eta = 1 - 1.5 eps
RHO="${RHO:-0.5}"
MC_CORES="${MC_CORES:-20}"
ARRAY_LIMIT="${ARRAY_LIMIT:-10}"

mkdir -p "$OUTDIR"

echo "----------------------------------------"
echo "Submitting robustness simulation"
echo "MODE:        $MODE"
echo "GRID:        $GRID"
echo "OUTDIR:      $OUTDIR"
echo "Tasks:       $NCELLS"
echo "R_REPS:      $R_REPS"
echo "R_BLOCK:     $R_BLOCK"
echo "B_PERM:      $B_PERM"
echo "DELTA_ROB:   $DELTA_ROB"
echo "ETA rule:    adaptive eta = 1 - 1.5 eps"
echo "ARRAY_LIMIT: $ARRAY_LIMIT"
echo "----------------------------------------"

START=1
while [ $START -le $NCELLS ]; do
  END=$((START + MAX_ARRAY - 1))
  if [ $END -gt $NCELLS ]; then
    END=$NCELLS
  fi

  echo "Submitting array chunk: ${START}-${END}"

  sbatch --array=${START}-${END}%${ARRAY_LIMIT} \
    --export=ALL,GRID_FILE="$GRID",STUDY="robustness",OUTDIR="$OUTDIR",BASE_SEED="$BASE_SEED",R_REPS="$R_REPS",R_BLOCK="$R_BLOCK",B_PERM="$B_PERM",ALPHA_LEVEL="$ALPHA_LEVEL",ETA="$ETA",DELTA_ROB="$DELTA_ROB",RHO="$RHO",MC_CORES="$MC_CORES" \
    scripts/slurm_hdmanova_array.sbatch

  START=$((END + 1))
done
