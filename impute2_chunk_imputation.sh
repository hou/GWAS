#!/bin/bash

if [ $# != 3 ];then
    echo "Usage: impute2_chunk_imputation.sh chunksize chr gwas_file_root_name"
    exit
fi

#Reference based chunking 
chunksize=$1
overlap=250000
chr=$2
gwas=$3
chrombegin=`zcat /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz | head -n 2 | tail -n 1 | cut -d" " -f2`
chromend=`zcat /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz | tail -n 1 | cut -d" " -f2`

for start in `seq $chrombegin $chunksize $chromend`
    do
        end=`expr $start + $chunksize - 1`
        dist2chromend=`expr $chromend - $start`
        if [ $dist2chromend -lt 7000000 ]; then
            echo "impute2 -use_prephased_g -m /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt -h /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz -l /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz -known_haps_g ${gwas}_chr$chr.phased.haps -o_gz -o ${gwas}_chr${chr}_impute2.$start.$chromend -int $start $chromend -filt_rules_l 'EUR==0'"
            break
        else 
            echo "impute2 -use_prephased_g -m /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt -h /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz -l /data/houl3/Share/Imputation/1000G/Impute2_reference/1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz -known_haps_g ${gwas}_chr$chr.phased.haps -o_gz -o ${gwas}_chr${chr}_impute2.$start.$end -int $start $end -filt_rules_l 'EUR==0'"
        fi
    done

#for chr in {1..22}; do bash impute2_chunk_imputation.sh 5000000 $chr WTCCC1_Cleaned; done >impute2.swarm
