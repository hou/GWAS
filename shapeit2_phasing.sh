#!/bin/bash

if [ ! $# -eq 1 ]; then
    echo "Usage: shapeit2_phasing.sh <root name of the GWAS files>"
    exit
fi

file=$1
for chr in {1..22}
do
    echo "shapeit2 -B ${file}_hg19_for_phasing_fwd_chr$chr -M /data/houl3/Share/Genetic_map/chr$chr.gmap.gz -O ${file}_chr$chr.phased -L ${file}_chr${chr}_shapeit2_phasing.log -T 16"
done
