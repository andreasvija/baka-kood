#!/bin/bash
#screen, ssh stage1

#SBATCH -p main
#SBATCH -J run_qtlmap_annots
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=24G

module load java-1.8.0_40
module load singularity
module load nextflow

git clone https://github.com/kerimoff/qtlmap.git

/gpfs/hpc/home/andreasv/baka/nextflow/nextflow \
	run /gpfs/hpc/home/andreasv/baka/pipeline/qtlmap/main_multi_study.nf \
	-profile tartu_hpc \
	--studyFile /gpfs/hpc/home/andreasv/baka/analysis/sources_annots.tsv \
	--cis_window 200000 \
	--run_permutation true \
	--is_imputed FALSE \
	-resume
