#!/bin/bash
#SBATCH --partition=cpu-preempt
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=2
#SBATCH --mem=2000
#SBATCH --time=15

#SBATCH --array=0 #1-204 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=../logs/log_2025_05_16_intersect_nlcd_mean/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=intersectNLCD

#setup env
module load conda/latest
conda activate narrowPIXC

#run code
python3 intersectNLCD.py mean