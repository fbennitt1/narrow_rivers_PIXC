#!/bin/bash
#SBATCH --partition=cpu-preempt # ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=10000 #64512  #
#SBATCH --time=20 #4320 #2880 #in minutes

#SBATCH --array=57,61 # 0-204 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=../logs/log_2025_04_10_buffer_flow_mean/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=bufferFlowlines

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
dates

#run code
python3 bufferFlowlines.py mean

#time stamp
date