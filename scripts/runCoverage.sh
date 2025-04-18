#!/bin/bash
#SBATCH --partition=cpu-preempt #ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20000
#SBATCH --time=15 #in minutes 

#SBATCH --array=3203 #3109-4662 # 0-1554 # 1555-3108
#SBATCH --output=./log_2025_03_25_HUC2_01_mean_profile/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=evalCoverage

#setup env
module load conda/latest
conda activate narrowPIXC

#run code 
python3 evalCoverage.py mean # CHANGE