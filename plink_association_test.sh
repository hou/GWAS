#!/bin/bash

if [ ! $# -eq 2 ]; then
    echo "Usage: plink_association_test.sh <root name of the GWAS files> <No. of chunks>"
    exit
fi

file=$1
n_chunks=$2

if [ -e "${file}_INFO_ge0.5.map" ]
    then
    rm -rf ${file}_INFO_ge0.5.map
    echo "Found file '${file}_INFO_ge0.5.map', need to remove it first!"
fi

echo "Making map file based on impute2 info files..."
echo "Will only keep markers with a INFO >= 0.5 & MAF >= 0.01"
for chr in {1..22}
do
chunks=`ls ${file}_chr${chr}_impute2.*_info | sort -V`
for i in $chunks
do
awk -v chrom=$chr '{if (NR>1 && $5>=0.5 && $4>=0.01 && $4<=0.99) print chrom, $2, 0, $3}' $i >>${file}_INFO_ge0.5.map
done
echo "Chromosome $chr: Done!"
done
n=`wc -l ${file}_INFO_ge0.5.map | cut -d" " -f1`
echo "Done, there are $n markers in the map file!"
n_part=`echo "$n/$n_chunks + 1" | bc`
split -l $n_part -d ${file}_INFO_ge0.5.map

echo "Making dosage file list for PLINK..."
ls ${file}_chr*_info | sed -e 's/_info/.gz/' | sort -V >impute2_prob.list
awk '{if (NR>2) print $1,$2}' ${file}_chr20.phased.sample >${file}_impute2.ind
awk '{if (NR>2) print $1,$2,$4,$5,$6,$7}' ${file}_chr20.phased.sample >${file}_final.fam
awk -v prefix=$file '{print NR,$1,prefix"_impute2.ind"}' impute2_prob.list >plink_dosage.list
m=`wc -l plink_dosage.list | cut -d" " -f1`
echo "Done, there are $m dosage files in the list file."

if [ -e "plink_commands.swarm" ]
    then
    rm -rf plink_commands.swarm
    echo "Found file 'plink_commands.swarm', need to remove it first!"
fi

for i in {0..9}
do
    echo "plink --fam ${file}_final.fam --map x0$i --dosage plink_dosage.list list sepheader Zin format=3 skip0=1 skip1=1 --covar ${file}_PCA.cov --out ${file}_1000g_imputed_part0$i --noweb" >>plink_commands.swarm
done

end=`expr $n_chunks - 1`
for i in $(eval echo "{10..$end}")
do
    echo "plink --fam ${file}_final.fam --map x$i --dosage plink_dosage.list list sepheader Zin format=3 skip0=1 skip1=1 --covar ${file}_PCA.cov --out ${file}_1000g_imputed_part$i --noweb" >>plink_commands.swarm
done

#echo "plink --fam ${file}_final.fam --map ${file}_INFO_ge0.5.map --dosage plink_dosage.list list sepheader Zin format=3 skip0=1 skip1=1 --out ${file}_1000g_imputed --noweb" >plink.swarm
echo "The 'plink_commands.swarm' is ready, please submit it to biowulf."
echo "Good luck, bye..."
