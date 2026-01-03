#!/bin/bash
rm -rf /home/bt6725/tregs/logs_opt
rm -rf /home/bt6725/tregs/spsa_*
rm -rf /home/bt6725/tregs/merged.txt
sbatch /home/bt6725/tregs/DLL_run_opt_sub.sh
