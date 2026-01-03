#!/bin/bash
rm -rf /home/bt6725/tregs/logs_abm
rm -rf /scratch/gpfs/CMETCALF/sim_abm
sbatch /home/bt6725/tregs/DLL_run_abm_sub.sh
