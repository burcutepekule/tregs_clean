#!/bin/bash
rm -rf /home/bt6725/tregs/logs_opt
rm -rf /home/bt6725/tregs/spsa_*
rm -rf /home/bt6725/tregs/random_*
rm -rf /home/bt6725/tregs/merged_*
# sbatch /home/bt6725/tregs/DLL_run_opt_sub_1300.sh
# sbatch /home/bt6725/tregs/DLL_run_opt_sub_60800.sh
sbatch /home/bt6725/tregs/DLL_run_opt_sub_5504.sh

