#submit batch jobs in align
#use screen
#ssh stage1
#module load python-3.6.0
#python script requires a created SlurmOut folder

#dry run
#snakemake -s all_bwa.Snakefile -n -p out.txt --configfile config.yaml > drysnakeoutput
#visualize
#snakemake -s all_bwa.Snakefile --dag out.txt --configfile config.yaml | dot -Tsvg > dag.svg
#run in cluster
#snakemake -s all_bwa.Snakefile --cluster ../snakemake_submit_UT.py -p out.txt --configfile config.yaml --jobs 8

rule make_all:
	input:
		expand("results/{sample}.sorted.bam", sample=config["samples"]),
		expand("results/{sample}.sorted.bam.bai", sample=config["samples"])
	output:
		"out.txt"
	resources:
		mem = 1000
	threads: 1
	shell:
		"echo 'Done!' >> {output}"

rule align_reads:
	input:
		inputfile = "/gpfs/hpc/projects/genomic_references/Garieri_2017/{sample}.fastq.gz"
	output:
		outputfile = "results/{sample}.bam"
	params:
		index = "/gpfs/hpc/projects/genomic_references/annotations/GRCh38/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa",

		rg = "@RG\tID:{sample}\tSM:{sample}",

		tempfq = "temps/{sample}tempfq.fastq.gz",
		tempsai = "temps/{sample}tempsai.sai",
		tempbam = "temps/{sample}tempbam.bam"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load bwa-0.7.12
		module load samtools-1.9
		
		cp {input.inputfile} {params.tempfq}

		bwa aln -t {threads} {params.index} {params.tempfq} > {params.tempsai}
		bwa samse -r '{params.rg}' {params.index} {params.tempsai} {params.tempfq} | samtools view -b - > {params.tempbam}

		cp {params.tempbam} {output.outputfile}
		
		rm {params.tempfq}
		rm {params.tempsai}
		rm {params.tempbam}
		
		echo 'Completed aligning {wildcards.sample}' >> progress.txt
		"""

rule sort_and_index:
	input:
		inputfile = "results/{sample}.bam"
	output:
		outputfile = "results/{sample}.sorted.bam",
		outputindex = "results/{sample}.sorted.bam.bai"
	params:
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load samtools-1.9

		samtools sort -o {output.outputfile} {input.inputfile}
		samtools index {output.outputfile} {output.outputindex}

		rm {input.inputfile}
		
		echo 'Completed indexing {wildcards.sample}' >> progress.txt
		"""
