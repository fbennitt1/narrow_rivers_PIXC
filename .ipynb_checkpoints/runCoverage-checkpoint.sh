#!/bin/bash
#SBATCH --partition=cpu-preempt #ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20000
#SBATCH --time=15 #in minutes 

#SBATCH --array=3109-4662 # 0-1554 # 1555-3108 #
#SBATCH --output=./log_2025_03_02_max/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=evalCoverage

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
date

#run code 
python3 evalCoverage.py max # CHANGE

#time stamp
date