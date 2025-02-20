#!/bin/bash
#SBATCH --partition=cpu-preempt #ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=15000
#SBATCH --time=10 #in minutes 

#SBATCH --array=3203 #3981-5963 #1999-3980 #0-1998
#SBATCH --output=./log_2025_02_14_profiling_max/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=evalCoverage

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
date

#run code 
python3 evalCoverage.py max

#time stamp
date