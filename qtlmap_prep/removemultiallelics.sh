#!/bin/bash

#SBATCH -p main
#SBATCH -J splitmultiallelics
#SBATCH -t 07:00:00
#SBATCH -c 4
#SBATCH --mem=8G

module load bcftools-1.9

bcftools view --max-alleles 2 /gpfs/hpc/home/a72094/datasets/controlled_access/Garieri_2017/Garieri_filtered.vcf.gz > Garieri_filtered_2.vcf 
gzip Garieri_filtered_2.vcf
