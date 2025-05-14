#!/bin/bash
#SBATCH --partition=cpu-preempt # ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000 #64512  #
#SBATCH --time=20 #4320 #2880 #in minutes

#SBATCH --array=197 #0-204 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=../logs/log_2025_05_14_buffer_flow_nlcd_max/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=bufferFlowlinesNLCD

#setup env
module load conda/latest
conda activate narrowPIXC

#time stamp
date

#run code
python3 bufferFlowlinesNLCD.py max

#time stamp
date