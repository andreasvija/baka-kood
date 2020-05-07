#!/bin/bash

#SBATCH -p main
#SBATCH -J one_bam_to_bigwig
#SBATCH -t 07:00:00
#SBATCH -c 4
#SBATCH --mem=16G

sample="ERR2021349"
infile="results/${sample}.sorted.bam"
intermediate="temps/${sample}.bedgraph"
outfile="results/${sample}.bw"

module load bedtools/2.27.0
bedtools genomecov -bg -ibam $infile > $intermediate
C_COLLATE=C sort -k1,1 -k2,2n $intermediate -o $intermediate

/gpfs/hpc/home/andreasv/bedGraphToBigWig/bedGraphToBigWig $intermediate hg38.chrom.sizes $outfile
rm $intermediate
