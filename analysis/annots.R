wd = "/gpfs/hpc/home/andreasv/baka/analysis/"
#wd = "/home/kasutaja/Desktop/baka/analysis/"
#wd = "C:/Proge/Ülikool/baka/analysis/"
setwd(wd)

library("data.table") # %like%
library("rtracklayer")
library("GenomicRanges")
library("GenomicFeatures")
library("wiggleplotr")

library("dplyr")
library("ggplot2")
library("readr")
library("stringr")
library("reshape2")
library("tidyr")

PROMOTER_BUFFER_RANGE = 25
EXON_START_BUFFER_RANGE = 25
EXON_BACK_RANGE = 100


col_names = c("gene", "chr", "phenotype_start", "phenotype_end", "strand", "top_phenotype", "total_phenotypes", 
              "variants_tested", "variant_distance", "top_variant", "top_variant_chr", "top_variant_start", 
              "top_variant_end", "df", "dummy", "beta_1", "beta_2", "p_nominal", "slope", "p_adjusted", "p_adjusted_beta")
col_types = "cciicciiicciiiddddddd"

bwa_results = read_delim("bwa.permuted.txt", delim=" ", col_types=col_types, col_names=col_names)

txrevise_results = read_delim("txrev.permuted.txt", delim=" ", col_types=col_types, col_names=col_names)
fun = function(x, n) {
  return(strsplit(x, ".", fixed=TRUE)[[1]][n])
}
txrevise_results = txrevise_results %>%
  mutate(kind=sapply(gene, fun, n=2), gene=sapply(gene, fun, n=1)) %>%
  filter(kind=="upstream")


shared_genes = intersect(bwa_results$gene, txrevise_results$gene)
shared_genes = shared_genes

promoter_annots = read_tsv("../qtlmap_prep/FANTOM5_promoter_annotations.tsv", col_types="ciciciiccii") %>%
  filter(gene_name %in% shared_genes) %>%
  select(gene_name, peak_start, peak_end, strand, tss_id, chr)


upstream1 = GenomicFeatures::makeTxDbFromGFF("txrevise.grp_1.upstream.gff3")
upstream2 = GenomicFeatures::makeTxDbFromGFF("txrevise.grp_2.upstream.gff3")
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(upstream2, by = "tx", use.names = TRUE)


#for every gene shared by cage and txrevise (N = 9 084)

all_new_transcripts = list()
new_transcript_genes = c()

start_time = Sys.time()
for (gene in shared_genes) { # c("ENSG00000151694")
  gene_new_transcripts = list()
  
  #for every promoter belonging to the gene (N = 53 364)
  
  promoters = promoter_annots %>%
    filter(gene_name == gene)
  
  utilized_promoters = list() # have been used for creating transcripts
  for (i in 1 : dim(promoters)[1]) {
    
    peak_start = promoters$peak_start[i]
    peak_end = promoters$peak_end[i]
    strand = promoters$strand[i]
    chr = promoters$chr[i]
    
    #if promoter +- PROMOTER_BUFFER_RANGE does not overlap with an already utilized promoter
    
    promoter = GRanges(seqnames=chr, IRanges(start=peak_start, end=peak_end))
    extended_promoter = GRanges(seqnames=chr, IRanges(start=peak_start - PROMOTER_BUFFER_RANGE, 
                                                    end=peak_end + PROMOTER_BUFFER_RANGE))
    
    intersections = lapply(utilized_promoters, pintersect, extended_promoter, drop.nohit.ranges=TRUE)
    if (sum(unlist(lapply(intersections, length))) != 0) {next}
    
    #for every exon in transcripts (N = 109 847)
    
    annot1 = exons_list1[names(exons_list1) %like% paste0(gene, ".grp_1.upstream")]
    annot2 = exons_list2[names(exons_list2) %like% paste0(gene, ".grp_2.upstream")]
    exons = c(annot1, annot2)
    
    overlaps_with_exon_start = FALSE
    for (j in 1 : length(exons@unlistData)) {
      
      #if promoter +- EXON_START_BUFFER_RANGE does not overlap with an exon's start
      a = peak_start - EXON_START_BUFFER_RANGE
      b = peak_end + EXON_START_BUFFER_RANGE
      
      exon_start = exons@unlistData@ranges@start[j]
      if (strand == "-") {
        exon_start = exon_start + exons@unlistData@ranges@width[j]
      }
      
      if (a < exon_start & b > exon_start) {
        overlaps_with_exon_start = TRUE
        break
      }
      
    }
    if (overlaps_with_exon_start) {next}
    
    #for every transcript belonging to the gene
    
    for (j in 1 : length(ranges(exons))) {
      transcript = unlist(ranges(exons)[j])
      
      #for every exon belonging to the transcript
      
      for (k in 1 : length(transcript)) {
        exon_start = transcript@start[k]
        exon_end = transcript@start[k] + transcript@width[k]
        
        #if promoter overlapping exon + EXON_BACK_RANGE bp back
        
        if (strand == "+") {exon_start = exon_start - EXON_BACK_RANGE}
        if (strand == "-") {exon_end = exon_end + EXON_BACK_RANGE}
        if ((peak_start > exon_end) | (peak_end < exon_start)) {next}
        
        #create a new exon from start of promoter to end of overlapping exon
        
        if (strand == "+") {
          new_exon_start = peak_start
          new_exon_end = exon_end
        }
        if (strand == "-") {
          new_exon_start = exon_start
          new_exon_end = peak_end
        }
        
        new_exon = IRanges(start=new_exon_start, end=new_exon_end)
        
        #take all exons in the transcript after the overlapping exon and create new transcript
        
        following = IRanges()
        
        if (strand == "+") {
          if (k < length(transcript)) {
            following = unname(transcript[(k+1) : length(transcript)])
          }
          new_transcript = c(new_exon, following)
        }
        
        if (strand == "-") {
          if (k > 1) {
            following = unname(transcript[1 : (k-1)])
          }
          new_transcript = c(following, new_exon)
        }
        
        ns = 1 : length(new_transcript)
        exon_rank = ns
        if (strand == "-") {exon_rank = rev(exon_rank)}
        
        new_transcript = GRanges(seqnames=chr, ranges=new_transcript, strand=strand, 
                                 exon_id=ns, exon_name=ns, exon_rank=exon_rank) # Parent, gene_id ?
        new_is_duplicate = FALSE
        if (length(gene_new_transcripts) > 0) {
          new_is_duplicate = any(sapply(new_transcript %in% gene_new_transcripts, all))
        }
        
        
        if (!new_is_duplicate) {
          gene_new_transcripts = c(gene_new_transcripts, list(new_transcript))
          new_transcript_name = paste(gene, "new", "upstream", length(all_new_transcripts)+1, sep=".")
          all_new_transcripts[new_transcript_name] = new_transcript
          new_transcript_genes = c(new_transcript_genes, gene)
        }
        
      }
    }
    
    utilized_promoters = c(utilized_promoters, list(promoter))
  }
  
  # visualization
  '
  spec = promoters
  
  spec_wide = spec %>%
    mutate(peak_start = peak_start - 150, peak_end = peak_end + 150)
  
  rangeslists = c()
  for (i in 1:dim(spec)[1]) {
    row = spec[i,]
    rangeslists = c(rangeslists, 
                    GRanges(seqnames=2, 
                            ranges=IRanges(as.numeric(spec[i,2]), as.numeric(spec[i,3])), 
                            strand=as.character(spec[i,4])))
  }
  names(rangeslists) = spec$tss_id
  
  rangeslists_wide = c()
  for (i in 1:dim(spec_wide)[1]) {
    row = spec_wide[i,]
    rangeslists_wide = c(rangeslists_wide, 
                         GRanges(seqnames=2, 
                                 ranges=IRanges(as.numeric(spec_wide[i,2]), as.numeric(spec_wide[i,3])), 
                                 strand=as.character(spec_wide[i,4])))
  }
  names(rangeslists_wide) = spec_wide$tss_id
  
  
  filter_start = min(spec$peak_start) - 100000
  filter_end = max(spec$peak_end) + 100000
  region_filter = GRanges(seqnames=2, IRanges(start=filter_start, end=filter_end))
  
  filtered_exons = lapply(exons, pintersect, region_filter, drop.nohit.ranges=TRUE)
  filtered_exons = filtered_exons[lapply(filtered_exons, length) > 0]
  #filtered_exons = pintersect(GRanges(seqnames=2, transcript), region_filter, drop.nohit.ranges=TRUE)
  
  all_trs = c(filtered_exons, rangeslists, gene_new_transcripts)
  expanded_trs = c(filtered_exons, rangeslists_wide, gene_new_transcripts)
  
  
  print(plotTranscripts(exons=expanded_trs, cdss=all_trs))
  '
}

Sys.time() - start_time
length(all_new_transcripts)

all_new_transcripts = GRangesList(all_new_transcripts)
export.gff3(all_new_transcripts, "new_transcripts_25.gff3")

saveRDS(all_new_transcripts, "new_transcripts_25.rds")
saveRDS(new_transcript_genes, "new_transcript_genes_25.rds")
