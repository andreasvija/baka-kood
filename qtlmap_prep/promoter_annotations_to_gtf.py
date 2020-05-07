"""
tss_id
promoter_rank
gene_name
gene_count
chr
peak_start
peak_end
peak_score
strand
tss_start
tss_end
->
seqname - chr
source - fantom5
feature - promotor
start - peak_start
end - peak_end
score - peak_score
strand - strand
frame - .
attribute:
gene_id gene_name; // meta-feature
promotor_id tss_id; // feature
"""

inputfilename = 'FANTOM5_promoter_annotations.tsv'
outputfilename = 'FANTOM5_promoter_annotations.gtf'

with open(inputfilename, 'r') as infile:
	with open(outputfilename, 'w') as outfile:

		infile.readline()
		for inline in infile:
			#print(inline)
			
			incols = inline.strip().split('\t')
			outline = ''
			
			outline += incols[4] + '\t' # seqname - chr
			outline += 'fantom5' + '\t' # source - fantom5
			outline += 'promotor' + '\t' # feature - promotor
			outline += incols[5] + '\t' # start - peak_start
			outline += incols[6] + '\t' # end - peak_end
			outline += incols[7] + '\t' # score - peak_score
			outline += incols[8] + '\t' # strand - strand
			outline += '.' + '\t' # frame - .
			
			outline += 'gene_id' + ' "' + incols[0] + '";' # promotor_id tss_id; // feature
			outline += 'geeni_id' + ' "' + incols[2] + '";\t' # gene_id gene_name; // meta-feature
			outline += '\n'
			
			#print(outline)
			outfile.write(outline)
