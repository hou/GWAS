#!/bin/bash

if [ ! $# -eq 2 ]; then
    echo "Usage: dataPrep4Phasing.sh <root name of the GWAS files> <genome_build>"
    echo "For example: dataPrep4Phasing.sh WTCCC_QCed hg17/18/19"
    exit
fi

file=$1
build=$2
echo "Started: `date`"
if [ $build != "hg19" ]
    then
        awk '{print "chr"$1,$4-1,$4,$2}' ${file}.bim >${file}_$build.bed
        sed -i -e 's/chr23/chrX/' ${file}_$build.bed 
        sed -i -e 's/chr24/chrY/' ${file}_$build.bed
        sed -i -e 's/chr25/chrX/' ${file}_$build.bed
        sed -i -e 's/chr26/chrM/' ${file}_$build.bed
        echo "liftOver from $build to hg19"
        liftOver ${file}_$build.bed /data/houl3/Share/liftOver/${build}ToHg19.over.chain ${file}_hg19.bed ${file}_hg19.unMapped
        awk '{print $4,$3}' ${file}_hg19.bed > ${file}_hg19.pos
        grep -e '#' -v ${file}_hg19.unMapped | awk '{print $4}' >${file}_hg19.unMapped.snps
        n=`wc -l ${file}_hg19.unMapped.snps | cut -d" " -f1`
        echo "$n SNPs in the ${file}_hg19.unMapped.snps file will be excluded after liftOver!!!"
        plink --bfile ${file} --exclude ${file}_hg19.unMapped.snps --update-map ${file}_hg19.pos --make-bed --out ${file}_liftTohg19 --noweb
        Rscript ~/bin/bimQC.R ${file}_liftTohg19.bim autosomes noATorCG ${file}_hg19_autosomes_noATorCG.bim
        plink --bfile ${file}_liftTohg19 --extract ${file}_hg19_autosomes_noATorCG.bim --make-bed --out ${file}_hg19_autosomes --noweb
    else
        echo "Yeah! the data is already at $build, liftOver is not needed :)"
        Rscript ~/bin/bimQC.R ${file}.bim autosomes noATorCG ${file}_hg19_autosomes_noATorCG.bim
        plink --bfile ${file} --extract ${file}_hg19_autosomes_noATorCG.bim --make-bed --out ${file}_hg19_autosomes --noweb
fi
pseq0.08 $file new-project --resources /data/houl3/Share/hg19/
pseq0.08 $file load-plink --check-reference --fix-strand --file ${file}_hg19_autosomes --id Bipolar
pseq0.08 $file write-ped --name  ${file}_hg19_for_phasing_fwd --family-id
plink --tfile ${file}_hg19_for_phasing_fwd --make-bed --out ${file}_hg19_for_phasing_fwd --noweb
mv ${file}_hg19_for_phasing_fwd.fam ${file}_hg19_for_phasing_fwd.fam.archive
cp ${file}_hg19_autosomes.fam ${file}_hg19_for_phasing_fwd.fam
for chr in {1..22}
do
    plink --bfile ${file}_hg19_for_phasing_fwd --chr $chr --make-bed --out ${file}_hg19_for_phasing_fwd_chr$chr --noweb
done
echo "Finished: `date`"
