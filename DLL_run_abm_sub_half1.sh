#!/bin/bash
#SBATCH --job-name=treg_array_1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --array=0-197
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_abm/treg_1_%A_%a.out
#SBATCH --error=logs_abm/treg_1_%A_%a.err

module load anaconda3/2023.3
conda activate env_Treg

N_CHUNKS=198

CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running chunk $CHUNK_ID on node $(hostname)"

Rscript /home/bt6725/tregs/DLL_datagen_abm_half1.R \
    $N_CHUNKS \
    $CHUNK_ID

echo "Chunk $CHUNK_ID completed"
