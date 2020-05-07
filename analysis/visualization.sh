#!/bin/bash

#SBATCH -p main
#SBATCH -J run_visualizations
#SBATCH -t 23:00:00
#SBATCH -c 4
#SBATCH --mem=16G

module load R/3.6.1
R < visualization.R --no-save
