#!/bin/bash
#submit batch jobs in align
#use screen
#ssh stage1
#module load python-3.6.0
#python script requires a created SlurmOut folder

#dry run
#snakemake -s bam_to_bigwig.Snakefile -n -p out_bigwig.txt --configfile config.yaml > drysnakeoutput
#visualize
#snakemake -s bam_to_bigwig.Snakefile --dag out_bigwig.txt --configfile config.yaml | dot -Tsvg > dag.svg
#run in cluster
#snakemake -s bam_to_bigwig.Snakefile --cluster ../snakemake_submit_UT.py -p out_bigwig.txt --configfile config.yaml --jobs 8

rule make_all:
	input:
		expand("results/{sample}.bw", sample=config["samples"]),
		expand("results/{sample_geuvadis}.bw", sample_geuvadis=config["samples_geuvadis"])
	output:
		"out_bigwig.txt"
	resources:
		mem = 1000
	threads: 1
	shell:
		"echo 'Done!' >> {output}"

rule convert_geuvadis:
	input:
		infile = "results/{sample}.sorted.bam"
	output:
		outfile = "results/{sample}.bw"
	params:
		intermediate="temps/{sample}.bedgraph",
		chrom_sizes="hg38.chrom.sizes"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load bedtools/2.27.0
		bedtools genomecov -bg -ibam {input.infile} > {params.intermediate}
		C_COLLATE=C sort -k1,1 -k2,2n {params.intermediate} -o {params.intermediate}
		
		/gpfs/hpc/home/andreasv/bedGraphToBigWig/bedGraphToBigWig {params.intermediate} {params.chrom_sizes} {output.outfile}
		rm {params.intermediate}
		"""

rule convert_cage:
	input:
		infile = "/gpfs/hpc/projects/genomic_references/GEUVADIS/bams/{sample_geuvadis}.sorted.bam"
	output:
		outfile = "results/{sample_geuvadis}.bw"
	params:
		intermediate="temps/{sample_geuvadis}.bedgraph",
		chrom_sizes="hg38.chrom.sizes"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load bedtools/2.27.0
		bedtools genomecov -bg -ibam {input.infile} > {params.intermediate}
		C_COLLATE=C sort -k1,1 -k2,2n {params.intermediate} -o {params.intermediate}
		
		/gpfs/hpc/home/andreasv/bedGraphToBigWig/bedGraphToBigWig {params.intermediate} {params.chrom_sizes} {output.outfile}
		rm {params.intermediate}
		"""
