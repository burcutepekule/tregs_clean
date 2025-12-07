#!/bin/bash
#SBATCH --job-name=treg_array
#SBATCH --array=0-999
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs/treg_%A_%a.out
#SBATCH --error=logs/treg_%A_%a.err

module load anaconda3/2023.3
conda activate env_Treg

N_CHUNKS=1000

# SLURM_ARRAY_TASK_ID gives 0..99
CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

echo "Running chunk $CHUNK_ID on node $(hostname)"

Rscript /home/bt6725/tregs/DLL_datagen_abm.R \
    $N_CHUNKS \
    $CHUNK_ID

echo "Chunk $CHUNK_ID completed"
