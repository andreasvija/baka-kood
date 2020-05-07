#!/bin/bash

#SBATCH -p main
#SBATCH -J get_variant_info
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=8G

gunzip -c /gpfs/hpc/projects/GENCORD/Garieri_2017/Garieri_filtered.no_DS.vcf.gz > temp.vcf
gunzip -c /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz > temp_geuvadis.vcf

while read VARIANT
do
  cat temp.vcf | grep $VARIANT >> variantinfo.vcf
  cat temp_geuvadis.vcf | grep $VARIANT >> variantinfo_geuvadis.vcf
done < variants.txt

rm temp.vcf
rm temp_geuvadis.vcf
