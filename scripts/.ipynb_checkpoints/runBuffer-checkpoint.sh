#!/bin/bash
#SBATCH --partition=ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=15000 #64512  #
#SBATCH --time=15 #4320 #2880 #in minutes

#SBATCH --array=#1-204 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=../logs/log_2025_04_01_min/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=bufferNHD

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
date

#run code
python3 bufferNHD.py True min

#time stamp
date