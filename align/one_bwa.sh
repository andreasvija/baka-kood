#!/bin/bash

#SBATCH -p main
#SBATCH -J test_bwa_align
#SBATCH -t 07:00:00
#SBATCH -c 4
#SBATCH --mem=8G

threads=4
workfolder="/gpfs/hpc/home/a72094/projects/CAGE_promoters/"
rg="@RG\tID:{sample}\tSM:{sample}"

outputfile="${workfolder}out.bam"
inputfile="/gpfs/hpc/projects/genomic_references/Garieri_2017/ERR2021349.fastq.gz"
indexfile="/gpfs/hpc/projects/genomic_references/annotations/GRCh38/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

tempfq="${workfolder}tempfq.fastq.gz"
tempsai="${workfolder}tempsai.sai"
tempbam="${workfolder}tempbam.bam"

module load samtools-1.9
module load bwa-0.7.12

cp $inputfile $tempfq

#index already exists
bwa aln -t $threads $indexfile $tempfq > $tempsai
bwa samse -r $rg $indexfile $tempsai $tempfq | samtools view -b - > $tempbam

cp $tempbam $outputfile
rm $tempfq
rm $tempsai
rm $tempbam