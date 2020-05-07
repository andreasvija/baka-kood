wd = "/gpfs/hpc/home/andreasv/baka/analysis/"
#wd = "/home/kasutaja/Desktop/baka/analysis/"
#wd = "C:/Proge/Ülikool/baka/analysis/"
setwd(wd)

library("GenomicRanges")
library("rtracklayer")
library("wiggleplotr")
library("cowplot")
library("data.table")
library("GenomicFeatures")

library("dplyr")
library("ggplot2")
library("readr")
library("stringr")
library("reshape2")
library("tidyr")

set.seed(123)




sample_mapping = read_tsv("../qtlmap_prep/sampleMetadata.tsv") %>%
  select(sample_id, genotype_id)
sample_mapping_ = read_tsv("GEUVADIS_EUR.tsv") %>%
  select(sample_id, genotype_id)
overlap = intersect(sample_mapping$genotype_id, sample_mapping_$genotype_id) # n=78


genotypes = read_tsv("variantinfo.vcf")[c(3,10:163)] %>%
  melt(id.vars=c("ID"))
names(genotypes) = c("variant", "genotype_id", "alleles")
genotypes = genotypes %>%
  merge(sample_mapping) %>%
  select(-genotype_id)

mapping = list("0", "1", "1", "2")
names(mapping) = c("0|0", "0|1", "1|0", "1|1")
genotypes$alleles = mapping[genotypes$alleles]

genotypes_ = read_tsv("variantinfo_geuvadis.vcf")[c(3,10:454)] %>%
  melt(id.vars=c("ID"))
names(genotypes_) = c("variant", "genotype_id", "alleles")
genotypes_ = genotypes_ %>%
  filter(genotype_id %in% sample_mapping_$genotype_id) %>%
  merge(sample_mapping_) %>%
  select(-genotype_id)

mapping_ = list("G_0", "G_1", "G_1", "G_2")
names(mapping_) = c("0|0", "0|1", "1|0", "1|1")
genotypes_$alleles = mapping_[genotypes_$alleles]

genotypes = rbind(genotypes, genotypes_)
genotypes$alleles = unlist(genotypes$alleles)


summ = read_tsv("../qtlmap_prep/bwa.counts.summary", col_types = paste0(c("c", rep("n", 154)), collapse=""))
renamer <- function(name) {
  name = str_replace(name, "../align/results/", "")
  name = str_replace(name, ".sorted.bam", "")
  return(name)
}

colnames(summ) = unname(sapply(colnames(summ), renamer))
summ = t(summ[1,2:length(summ)])
summ = data_frame(sample_id=rownames(summ), N=summ)

summ_ = read_tsv("geuvadis_counts.txt")
summ_ = colSums(summ_[,-1])
summ_ = data_frame(sample_id=labels(summ_), N=summ_) %>%
  filter(sample_id %in% sample_mapping_$sample_id) # let's only look at european individuals

summ = rbind(summ, summ_)

sample_data = data_frame(sample_id=summ$sample_id, scaling_factor=summ$N/1000000) %>%
  mutate(bigWig = paste0("/gpfs/hpc/home/andreasv/baka/align/results/", sample_id, ".bw"))


to_visualize = read_tsv("bwa_better.tsv")

promoterMetadata = read_tsv("../qtlmap_prep/FANTOM5_promoter_annotations.tsv", col_types="ciciciiccii") %>%
  filter(gene_name %in% to_visualize$gene)


upstream1 = GenomicFeatures::makeTxDbFromGFF("txrevise.grp_1.upstream.gff3")
upstream2 = GenomicFeatures::makeTxDbFromGFF("txrevise.grp_2.upstream.gff3")
exons_list1 = GenomicFeatures::exonsBy(upstream1, by = "tx", use.names = TRUE)
exons_list2 = GenomicFeatures::exonsBy(upstream2, by = "tx", use.names = TRUE)




for (n in c(1: (dim(to_visualize)[1]) )) { # c(1: (dim(to_visualize)[1]) )
temp = tryCatch({
  
  print(n)
  observed_gene = to_visualize$gene[n]
  observed_variant = to_visualize$top_variant[n]
  
  
  spec = promoterMetadata %>% 
    filter(gene_name==observed_gene)
  
  spec_wide = spec %>%
    mutate(peak_start = peak_start - 150, peak_end = peak_end + 150)
  
  
  rangeslists = c()
  for (i in 1:dim(spec)[1]) {
    row = spec[i,]
    rangeslists = c(rangeslists, 
                    GRanges(seqnames=2, 
                            ranges=IRanges(as.numeric(spec[i,6]), as.numeric(spec[i,7])), 
                            strand=as.character(spec[i,9])))
  }
  names(rangeslists) = spec$tss_id
  
  rangeslists_wide = c()
  for (i in 1:dim(spec_wide)[1]) {
    row = spec_wide[i,]
    rangeslists_wide = c(rangeslists_wide, 
                         GRanges(seqnames=2, 
                                 ranges=IRanges(as.numeric(spec_wide[i,6]), as.numeric(spec_wide[i,7])), 
                                 strand=as.character(spec_wide[i,9])))
  }
  names(rangeslists_wide) = spec_wide$tss_id

  
  annot1 = exons_list1[names(exons_list1) %like% paste0(observed_gene, ".grp_1.upstream")]
  annot2 = exons_list2[names(exons_list2) %like% paste0(observed_gene, ".grp_2.upstream")]
  exons = c(annot1, annot2)
  

  filter_start = min(spec$peak_start) - 1000
  filter_end = max(spec$peak_end) + 1000
  region_filter = GRanges(seqnames=2, IRanges(start=filter_start, end=filter_end))
  
  filtered_exons = lapply(exons, pintersect, region_filter, drop.nohit.ranges=TRUE)
  filtered_exons = filtered_exons[lapply(filtered_exons, length) > 0]

  
  all_trs = c(filtered_exons, rangeslists)
  expanded_trs = c(filtered_exons, rangeslists_wide)
  
  
  genotypes_sub = genotypes %>%
    filter(variant == observed_variant) %>%
    select(sample_id, alleles)
  
  sample_data_this = sample_data %>%
    merge(genotypes_sub) %>%
    mutate(track_id=alleles, colour_group=alleles) %>%
    select(-alleles)
  
  sample_data_this_cage = sample_data_this %>%
    filter(track_id %in% mapping)
  sample_data_this_geuvadis = sample_data_this %>%
    filter(track_id %in% mapping_)
  
  
  a = plotCoverage(exons=expanded_trs, cdss=all_trs, track_data=sample_data_this_cage, #transcript_annotations=annots, 
                   fill_palette = c("#a1dab4", "#41b6c4", "#225ea8", "#a1dab4", "#41b6c4", "#225ea8"), 
                   transcript_label=FALSE, plot_fraction=0.2, return_subplots_list=FALSE, 
                   rescale_introns=TRUE, new_intron_length=100, heights=c(0.5, 0.5))
  
  b = plotCoverage(exons=expanded_trs, cdss=all_trs, track_data=sample_data_this_geuvadis, #transcript_annotations=annots_, 
                   fill_palette = c("#a1dab4", "#41b6c4", "#225ea8", "#a1dab4", "#41b6c4", "#225ea8"), 
                   transcript_label=FALSE, plot_fraction=0.2, return_subplots_list=FALSE, 
                   rescale_introns=TRUE, new_intron_length=100, heights=c(0.5, 0.5))
  
  combo = plot_grid(a, b, ncol=1, nrow=2)
  
  ggsave(paste("plots/", n, observed_gene, ".pdf", sep="_"), combo)

}, error = function(e){print(e)})
}
