#!/bin/bash
#SBATCH --partition=cpu-preempt #ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=50000
#SBATCH --time=60 #in minutes 

#SBATCH --array=7992-9375
#SBATCH --output=./log_2024_12_04/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=prepCoverage

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
date

#run code 
python3 evalCoverage.py

#time stamp
date

echo "Script completed."