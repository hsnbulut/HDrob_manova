#!/bin/bash
set -euo pipefail

module purge
module load apps/R/4.5.2-gcc-14.1.0

cd "$(dirname "$0")/.."

mkdir -p logs
mkdir -p results/raw_chunks/type1
mkdir -p results/raw_chunks/power
mkdir -p results/raw_chunks/robustness
mkdir -p results/tables

Rscript code/simulation/create_simulation_grids.R

submit_in_chunks () {
  STUDY=$1
  GRID=$2
  OUTDIR=$3

  NCELLS=$(($(wc -l < "$GRID") - 1))
  MAX_ARRAY=99

  echo "----------------------------------------"
  echo "Submitting $STUDY"
  echo "GRID:   $GRID"
  echo "OUTDIR: $OUTDIR"
  echo "Tasks:  $NCELLS"
  echo "----------------------------------------"

  START=1
  while [ $START -le $NCELLS ]; do
    END=$((START + MAX_ARRAY - 1))
    if [ $END -gt $NCELLS ]; then
      END=$NCELLS
    fi

    echo "Submitting array chunk: ${START}-${END}"

    sbatch --array=${START}-${END}%20 \
      --export=ALL,GRID_FILE="$GRID",STUDY="$STUDY",OUTDIR="$OUTDIR",BASE_SEED=20260419,R_REPS=1000,R_BLOCK=100,B_PERM=199,ALPHA_LEVEL=0.05,ETA=0.75,DELTA_ROB=4.0,RHO=0.5,MC_CORES=20 \
      scripts/slurm_hdmanova_array.sbatch
    START=$((END + 1))
  done
}

submit_in_chunks "type1" "results/tables/grid_type1.csv" "results/raw_chunks/type1"
submit_in_chunks "power" "results/tables/grid_power.csv" "results/raw_chunks/power"
submit_in_chunks "robustness" "results/tables/grid_robustness.csv" "results/raw_chunks/robustness"