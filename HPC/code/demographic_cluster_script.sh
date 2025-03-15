#!/bin/bash
#PBS -lwalltime=00:05:00
#PBS -lselect=1:ncpus=1:mem=1gb


module load R
echo "R is about to run"
cd /rds/general/user/yh4724/home/
Rscript Demographic_cluster.R
echo "R has finished running"
