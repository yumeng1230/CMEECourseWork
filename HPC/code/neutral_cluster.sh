#!/bin/bash
#PBS -lwalltime=12:00:00
#PBS -lselect=1:ncpus=1:mem=4gb


module load R

echo "R is about to run"
cd /rds/general/user/yh4724/home/
Rscript yh4724_HPC_2024_neutral_cluster.R
echo "R has finished running"
