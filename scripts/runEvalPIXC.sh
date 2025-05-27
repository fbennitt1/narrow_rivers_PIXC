#!/bin/bash
#SBATCH --partition=cpu-preempt #ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80000
#SBATCH --time=40 #in minutes 

#SBATCH --array=5106,8617

#SBATCH --output=../logs/log_2025_05_27_HUC2_17_mean_1_redo/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=evalPIXC

#setup env
module load conda/latest
conda activate narrowPIXC

#run code 
python3 evalPIXC.py mean # CHANGE