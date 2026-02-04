#!/bin/bash
#SBATCH --job-name=treg_array
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --array=0-39
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_abm/treg_%x_%A_%a.out
#SBATCH --error=logs_abm/treg_%x_%A_%a.err

module load anaconda3/2023.3
conda activate env_Treg

N_SCENARIO_CHUNKS=40
SCENARIO_CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running param chunk $PARAM_CHUNK_ID, scenario chunk $SCENARIO_CHUNK_ID on node $(hostname)"

Rscript /home/bt6725/tregs/DLL_datagen_abm_chunks_combined.R \
    $N_SCENARIO_CHUNKS \
    $SCENARIO_CHUNK_ID

echo "Chunk completed"