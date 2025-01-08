#!/bin/bash
#SBATCH --partition=cpu-preempt

#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000
#SBATCH --time=60 #in minutes 

#SBATCH --array=3997-5993
#SBATCH --output=./log_2024_12_04_taylor/log_%A_%a.log # job id: %A_ 
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