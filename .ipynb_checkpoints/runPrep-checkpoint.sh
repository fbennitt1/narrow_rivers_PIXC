#!/bin/bash
#SBATCH --partition=ceewater_cjgleason-cpu #gpu-preempt #gpu #cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=2 
#SBATCH --mem=10000 #64512  #
#SBATCH --time=50 #4320 #2880 #in minutes

#SBATCH --array=0,2-204 # 51-57,108-112,118-124,145-146,169-171,195-196 #0-209 # %5 means limit the runs with 5 arrays at a time
#SBATCH --output=./log/log_%A_%a.log # job id: %A_ 
#SBATCH --job-name=prepNHD

#setup env
module load miniconda/22.11.1-1
conda activate small-rivers-1

#time stamp
date

#run code 
# python3 prepNHD.py
python3 prepNHDkeepWaterbody.py

#time stamp
date