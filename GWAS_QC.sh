#!/bin/bash

if [ ! $# -eq 2 ]; then
    echo "Usage: GWAS_QC.sh <root name of the GWAS files> genome_build"
    echo
    echo "If you have a gwas data set in plink binary format named as my_gwas.bed my_gwas.bim my_gwas.fam,"
    echo "and the base pair positons in the bim file are based on the genome build hg18."
    echo "You can just use the following command to QC your data:"
    echo
    echo "bash GWAS_QC.sh my_gwas hg18"
    exit
fi

if [ "$2" != "hg19" -a "$2" != "hg18" ]
    then
    echo
    echo "Error: the genomic build can only be either hg18 or hg19, please lift your data over first!"
    echo
    exit
fi

file=$1
qc_folder='./GWAS/Files'
nids=`wc -l $file.fam | awk '{print $1}'`
nmarkers=`wc -l $file.bim | awk '{print $1}'`
echo "Found $nids individuals, $nmarkers markers from your gwas data set!"

echo "Working on the sex checking ..."
echo "Please make sure your gwas data set includes some markers on chrX!!!"
plink --bfile $file --check-sex --out $file --noweb
grep -e 'PROBLEM' -e 'STATUS' $file.sexcheck >$file.sexprobs
awk 'BEGIN{OFS="\t"}{if(NR>1) print $1,$2}' $file.sexprobs >$file.sexprobs.ind
nsexprobs=`wc -l $file.sexprobs.ind | awk '{print $1}'`
echo "Found $nsexprobs individuals, for whom the reported sex in the PED file does not match the estimated sex!"
echo

echo "Estimating the missing rate ..."
plink --bfile $file --missing --out $file --noweb

echo "Working on the heterozygosity estimation ..."
plink --bfile $file --het --out $file --noweb

echo "Working on the LD-based pruning ..."

if [ ! -f $qc_folder/high-LD-regions-hg19.txt ]
    then
    git clone https://github.com/hou/GWAS.git
fi

if [ "$2" == "hg19" ]
    then
    plink --bfile $file --exclude $qc_folder/high-LD-regions-hg19.txt --range --indep-pairwise 50 5 0.2 --out $file --noweb
    else
        plink --bfile $file --exclude $qc_folder/high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out $file --noweb
fi

echo "Merging your GWAS data with the HapMap phase3 data set ..."
awk '$1>22' $file.bim >$file.nonAutosomes.snps
plink --bfile $file --extract $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim --make-bed --out $file.hapmap-snps --noweb
plink --bfile $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps --bmerge $file.hapmap-snps.bed $file.hapmap-snps.bim $file.hapmap-snps.fam --extract $file.prune.in --exclude $file.nonAutosomes.snps --make-bed --out $file.hapmap3r2.pruned --noweb
if [ -f $file.hapmap3r2.pruned.missnp ]; then
    plink --bfile $file.hapmap-snps --flip $file.hapmap3r2.pruned.missnp --make-bed --out $file.hapmap-snps-flipped --noweb
    rm -rf $file.hapmap3r2.pruned.missnp
    plink --bfile $file.hapmap-snps-flipped --bmerge $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --extract $file.prune.in --exclude $file.nonAutosomes.snps --make-bed --out $file.hapmap3r2.pruned --noweb
    if [ -f $file.hapmap3r2.pruned.missnp ]; then
        plink --bfile $file.hapmap-snps-flipped --exclude $file.hapmap3r2.pruned.missnp --make-bed --out $file.hapmap-snps-flipped --noweb
        cat $file.hapmap3r2.pruned.missnp $file.nonAutosomes.snps >$file.allExcluded.snps
        plink --bfile $file.hapmap-snps-flipped --bmerge $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim $qc_folder/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --extract $file.prune.in --exclude $file.allExcluded.snps --make-bed --out $file.hapmap3r2.pruned --noweb
        fi
fi

echo "Working on the principal component analysis ..."
echo
echo "Please be patient, it might take quite a while to finish"
smartpca.perl -i $file.hapmap3r2.pruned.bed -a $file.hapmap3r2.pruned.bim -b $file.hapmap3r2.pruned.fam -o $file.hapmap3r2.pca -p $file.hapmap3r2.plot -e $file.hapmap3r2.eval -l $file.hapmap3r2.pca.log -k 10 -t 10 -m 0

echo
echo "Working on the relationship checking ..."
plink --bfile $file --extract $file.prune.in --exclude $file.nonAutosomes.snps --maf 0.05 --genome --out $file --noweb
awk '$10>0.1' $file.genome >$file.related.ids
nrel=`wc -l $file.related.ids | awk '{print $1}'`
echo "Found $nrel pairs that have a PI_HAT > 0.1!"
echo "Warning: If your data set includes individuals with different ethnicities, you might have to run the relationship checking for each ethnicity seperately!"
echo
echo "Generating the 'missing rate' vs. the 'heterozygosity rate' scatter plot ..."
echo "imiss=read.table(\"$file.imiss\",h=T)" >$file.imiss.vs.het.R
echo "imiss\$logF_MISS=log10(imiss[,6])" >>$file.imiss.vs.het.R
echo "write.table(imiss[which(imiss[,6]>0.03),1:2],\"$file.imiss.gt3percent.ind\",row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')" >>$file.imiss.vs.het.R
echo "het=read.table(\"$file.het\",h=T)" >>$file.imiss.vs.het.R
echo "het\$meanHet=(het\$N.NM. - het\$O.HOM.)/het\$N.NM." >>$file.imiss.vs.het.R
echo "het.outliers <- het[het\$meanHet < mean(het\$meanHet)-(3*sd(het\$meanHet)) | het\$meanHet > mean(het\$meanHet)+(3*sd(het\$meanHet)),]" >>$file.imiss.vs.het.R
echo "write.table(het.outliers[,1:2],\"$file.het.outliers.ind\",row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')" >> $file.imiss.vs.het.R
echo "if(!require(geneplotter)){" >>$file.imiss.vs.het.R
echo "  source('http://bioconductor.org/biocLite.R')">>$file.imiss.vs.het.R
echo "  biocLite('geneplotter')}">>$file.imiss.vs.het.R
echo "colors  <- densCols(imiss\$logF_MISS,het\$meanHet)" >>$file.imiss.vs.het.R
echo "pdf(\"$file.imiss-vs-het.pdf\")" >>$file.imiss.vs.het.R
echo "plot(imiss\$logF_MISS,het\$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5),pch=20, xlab='Proportion of missing genotypes', ylab='Heterozygosity rate',axes=F)" >>$file.imiss.vs.het.R
echo "axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)" >>$file.imiss.vs.het.R
echo "axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))" >>$file.imiss.vs.het.R
echo "abline(h=mean(het\$meanHet)-(3*sd(het\$meanHet)),col='RED',lty=2)" >>$file.imiss.vs.het.R
echo "abline(h=mean(het\$meanHet)+(3*sd(het\$meanHet)),col='RED',lty=2)" >>$file.imiss.vs.het.R
echo "abline(v=-1.522879, col='RED', lty=2)" >>$file.imiss.vs.het.R

R --no-save <$file.imiss.vs.het.R
echo
echo "Generating the principal components scatter plots ..."
echo "data=read.table(\"$file.hapmap3r2.pca.evec\",head=F,skip=1)" >$file.pca.plot.R
echo "data\$color <- 'black'" >>$file.pca.plot.R
echo "data\$color[data\$V12 %in% c('5','4')] <- 'PURPLE'" >>$file.pca.plot.R
echo "data\$color[data\$V12 %in% '3'] <- 'RED'" >>$file.pca.plot.R
echo "data\$color[data\$V12 %in% '6'] <- 'GREEN'" >>$file.pca.plot.R
echo "index=which(data\$V12 %in% c('Case','Control'))" >>$file.pca.plot.R
echo "pdf(\"$file.pca.pdf\")" >>$file.pca.plot.R 
echo "plot(data\$V2,data\$V3,pch=20,xlab='PC1', ylab='PC2',col=data\$color)" >>$file.pca.plot.R
echo "points(data\$V2[index],data\$V3[index],pch=20,col='#00001110')" >>$file.pca.plot.R
echo "legend('topright',c('CEU','CHB+JPT','YRI',\"$file\"),col=c('RED','PURPLE','GREEN','BLACK'),pch=20)" >>$file.pca.plot.R
echo "plot(data[,2:4],,col=data\$color)">>$file.pca.plot.R
echo "dev.off()">>$file.pca.plot.R
echo "pca <- data[data\$color %in% 'black',]" >>$file.pca.plot.R
echo "pc1.outliers <- pca[pca\$V2 < (mean(data\$V2[data\$color %in% 'RED']) - 6*sd(data\$V2[data\$color %in% 'RED'])) | pca\$V2 > (mean(data\$V2[data\$color %in% 'RED']) + 6*sd(data\$V2[data\$color %in% 'RED'])),]" >>$file.pca.plot.R
echo "pc2.outliers <- pca[pca\$V3 < (mean(data\$V3[data\$color %in% 'RED']) - 6*sd(data\$V3[data\$color %in% 'RED'])) | pca\$V3 > (mean(data\$V3[data\$color %in% 'RED']) + 6*sd(data\$V3[data\$color %in% 'RED'])),]" >>$file.pca.plot.R
echo "pca.outliers <- unique(rbind(pc1.outliers, pc2.outliers))" >>$file.pca.plot.R
echo "pca.outliers\$ID <- gsub(':','\t',pca.outliers\$V1)" >>$file.pca.plot.R
echo "write.table(pca.outliers\$ID,\"$file.pca.outliers.ind\",row.name=FALSE, col.names=F, quote=F, sep='\t')" >>$file.pca.plot.R

R --no-save <$file.pca.plot.R

echo "imiss <- read.table(\"$file.imiss\",head=T,as.is=T)" >$file.relcheck.R
echo "ids <- read.table(\"$file.related.ids\",head=T,as.is=T)" >>$file.relcheck.R
echo "imiss <- imiss[,c('FID','IID','F_MISS')]" >>$file.relcheck.R
echo "names(imiss) <- c('FID2','IID2','F_MISS2')" >>$file.relcheck.R
echo "ids <- merge(ids,imiss,sort=F)" >>$file.relcheck.R
echo "names(imiss) <- c('FID1','IID1','F_MISS1')" >>$file.relcheck.R 
echo "ids <- merge(ids,imiss,sort=F)" >>$file.relcheck.R
echo "ids\$FID <- ifelse(ids\$F_MISS1 > ids\$F_MISS2, ids\$FID1, ids\$FID2)" >>$file.relcheck.R 
echo "ids\$IID <- ifelse(ids\$F_MISS1 > ids\$F_MISS2, ids\$IID1, ids\$IID2)" >>$file.relcheck.R
echo "ids <- ids[,c('FID','IID','FID1','IID1','F_MISS1','FID2','IID2','F_MISS2','RT','EZ','Z0','Z1','Z2','PI_HAT')]" >>$file.relcheck.R
echo "write.table(ids,\"$file.related.ind.info\",row.names=F,quote=F,sep='\t')" >>$file.relcheck.R

R --no-save <$file.relcheck.R

cat $file.het.outliers.ind  $file.imiss.gt3percent.ind  $file.pca.outliers.ind $file.sexprobs.ind >$file.het.imiss.pca.sex.probs.ind
grep -v -f $file.het.imiss.pca.sex.probs.ind -w $file.related.ind.info | awk 'BEGIN{OFS="\t"}{if(NR>1) print $1,$2}' >$file.relcheck.probs.ind
cat $file.het.imiss.pca.sex.probs.ind $file.relcheck.probs.ind >$file.to.be.deleted.ind
ndel=`wc -l $file.to.be.deleted.ind | awk '{print $1}'`
echo
echo "Done! $ndel individuals need to be excluded from your gwas data!"
echo
echo "Please use the following plink command to get the QCed data set:"
echo "plink --bfile $file --exclude $file.to.be.deleted.ind --geno 0.02 --maf 0.01 --hwe 0.000001 --make-bed --out ${file}_QCed --noweb"

