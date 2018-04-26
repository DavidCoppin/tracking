#!/bin/bash
#SBATCH --job-name=clusters
#SBATCH --account=uoa00510
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --constraint=avx
#SBATCH --partition=debug

## load dependencies
module load OpenCV/3.4.0-gimkl-2017a-Python-2.7.13

## Test Madagascar
srun python testCmorph.py -d1 1999-01-01 -d2 1999-01-04 -lons 250:750 -lats 0:300 -suffix madag
