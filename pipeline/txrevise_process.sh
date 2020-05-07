#!/bin/bash
#screen, ssh stage1

#SBATCH -p main
#SBATCH -J run_txrevise_process
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=24G

module load java-1.8.0_40
module load singularity
module load nextflow

git clone https://github.com/eQTL-Catalogue/rnaseq
git clone https://github.com/kerimoff/qcnorm

NXF_VER=18.10.1 nextflow \
	run rnaseq/main.nf \
	-profile eqtl_catalogue \
	--readPathsFile /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/readPaths/GEUVADIS_readPaths.tsv \
	--unstranded \
	--skip_qc \
	--skip_multiqc \
	--skip_stringtie \
	--saveAlignedIntermediates \
	--run_txrevise \
	--txrevise_gffs '/gpfs/hpc/projects/genomic_references/annotations/txrevise/Homo_sapiens.GRCh38.96_CAGE_25bp/*.gff3' \
	-resume \
	-executor.queueSize 100

nextflow \
	run qcnorm/normalisation.nf \
	-profile tartu_hpc \
	-resume \
	--study_name new_annots \
	--quant_results_path /gpfs/hpc/home/andreasv/baka/pipeline/results \
	--sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GEUVADIS_EUR.tsv \
	--txrev_pheno_meta_path /gpfs/hpc/projects/genomic_references/annotations/txrevise/Homo_sapiens.GRCh38.96_CAGE_25bp/txrevise_Ensembl_96_CAGE_25bp_phenotype_metadata.tsv.gz \
	--skip_exon_norm \
	--skip_tx_norm \
	--skip_leafcutter_norm \
	--outdir new_annots
