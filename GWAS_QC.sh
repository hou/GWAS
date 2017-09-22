#!/bin/bash

# Author: Liping Hou
# Date: April 15, 2015
# Updated: Sept 20, 2017
# This is a simple script that can QC your GWAS data
# Input: Genotype data in plink binary format and hapmap phase3 data needed for PCA analysis
#        The hapmap phase 3 data could be downloaded from this link: https://github.com/hou/GWAS/tree/master/Files
# This script will do:
#   1. Check the gender 
#   2. Estimate missing rate at both the variant and individual level
#   3. Identify population outliers (non-Caucasian) by running PCA analysis
#   4. Identify related pairs and extract the one with higher missing rate for each pair


if [ ! $# -eq 2 ]; then
    echo "Usage: GWAS_QC.sh [root name of the GWAS files] [hg18/hg19_flag]"
    echo "Example: GWAS_QC.sh NIMH hg19"
    exit
fi

# Need to load eigensoft and plink (v1.9) first
module load eigensoft plink

file=$1
build=$2
qc_folder='/data/houl3/Share/GWAS_QC'

awk '$1>22' $file.bim >$file.nonAutosomes.snps
awk '$1==25' $file.bim >$file.PAR.snps
# Gender check
# If the dataset already contains an XY region run --check-sex otherwise run --split-x first then --check-sex
if [ -s "$file.PAR.snps" ]; then
    plink --bfile $file --check-sex --out $file
else
    plink --bfile $file --split-x $build no-fail --make-bed --out $file.splitX
    plink --bfile $file.splitX --check-sex --out $file
fi
grep -e 'PROBLEM' -e 'STATUS' $file.sexcheck >$file.sexprobs
awk 'BEGIN{OFS="\t"}{if(NR>1) print $1,$2}' $file.sexprobs >$file.sexprobs.ind
# Missingness
plink --bfile $file --missing --out $file
# Inbreeding coefficient estimate
plink --bfile $file --het --out $file
# LD-pruning
if [ build == "hg19" ]; then
    plink --bfile $file --exclude range $qc_folder/high-LD-regions-hg19.txt --indep-pairwise 50 5 0.2 --out $file
else
    plink --bfile $file --exclude range $qc_folder/high-LD-regions.txt --indep-pairwise 50 5 0.2 --out $file
fi
# Merge your genotype data with hapmap data
plink --bfile $file --extract $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt --make-bed --out ${file}.hapmap-snps
plink --bfile $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps --bmerge $file.hapmap-snps.bed $file.hapmap-snps.bim $file.hapmap-snps.fam --extract $file.prune.in --exclude $file.nonAutosomes.snps --make-bed --out $file.hapmap3r2.pruned
if [ -f $file.hapmap3r2.pruned-merge.missnp ]; then
    plink --bfile $file.hapmap-snps --flip $file.hapmap3r2.pruned-merge.missnp --make-bed --out $file.hapmap-snps-flipped
    rm -rf $file.hapmap3r2.pruned-merge.missnp
    plink --bfile $file.hapmap-snps-flipped --bmerge $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --extract $file.prune.in --exclude $file.nonAutosomes.snps --make-bed --out $file.hapmap3r2.pruned
    if [ -f $file.hapmap3r2.pruned-merge.missnp ]; then
        plink --bfile $file.hapmap-snps-flipped --exclude $file.hapmap3r2.pruned-merge.missnp --make-bed --out $file.hapmap-snps-flipped
        cat $file.hapmap3r2.pruned-merge.missnp $file.nonAutosomes.snps >$file.allExcluded.snps
        plink --bfile $file.hapmap-snps-flipped --bmerge $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --extract $file.prune.in --exclude $file.allExcluded.snps --make-bed --out $file.hapmap3r2.pruned
    fi
fi
# PCA with EIGENSOFT
smartpca.perl -i $file.hapmap3r2.pruned.bed -a $file.hapmap3r2.pruned.bim -b $file.hapmap3r2.pruned.fam -o $file.hapmap3r2.pca -p $file.hapmap3r2.plot -e $file.hapmap3r2.eval -l $file.hapmap3r2.pca.log -k 10 -t 10 -m 0
# Relationship checking
plink --bfile $file --extract $file.prune.in --exclude $file.nonAutosomes.snps --maf 0.05 --genome --out $file
awk '$10>0.1' $file.genome >$file.related.ids

# Plotting and extract related individuals with higher missing rate
Rscript ~/bin/GWAS_QC_plot.R $file

cat $file.het.outliers.ind  $file.imiss.gt3percent.ind  $file.pca.outliers.ind $file.sexprobs.ind >het.imiss.pca.sex.probs.ind
grep -v -f het.imiss.pca.sex.probs.ind -w $file.related.ind.info | awk 'BEGIN{OFS="\t"}{if(NR>1) print $1,$2}' >relcheck.probs.ind
cat het.imiss.pca.sex.probs.ind relcheck.probs.ind >to.be.deleted.ind

#plink --bfile $file --remove to.be.deleted.ind --maf 0.01 --geno 0.05 --hwe 0.0001 --make-bed --out $file.QCed
