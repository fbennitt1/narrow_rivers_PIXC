#!/bin/bash
#SBATCH --partition=cpu-preempt # ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=40000 #64512  #
#SBATCH --time=20 #4320 #2880 #in minutes

#SBATCH --array=197,108# 3,19,20,22,32,60 # 0-204 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=../logs/log_2025_04_30_buffer_mean/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=bufferNHD

#setup env
module load conda/latest
conda activate narrowPIXC

#run code
python3 bufferSegments.py max