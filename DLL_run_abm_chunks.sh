#!/bin/bash
rm -rf /home/bt6725/tregs/logs_abm
rm -rf /scratch/gpfs/CMETCALF/sim_abm

for i in {1..16}; do
    sbatch --job-name=treg_array_${i} /home/bt6725/tregs/DLL_run_abm_sub_chunk.sh $i
done