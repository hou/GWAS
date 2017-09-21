args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
    stop("Usage: Rscript GWAS_QC_plot.R <root name of the GWAS files>")
}else{
    options(echo=TRUE)
}

fileName = args[1]

# Genotype missing rate versus heterozygosity 
imiss.file = paste(fileName,".imiss",sep="")
imiss=read.table(imiss.file,head=TRUE)
imiss$logF_MISS = log10(imiss[,6])
imiss.output = paste(fileName,".imiss.gt3percent.ind",sep="")
write.table(imiss[which(imiss[,6]>0.03),1:2],imiss.output,row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
het.file = paste(fileName,".het",sep="")
het=read.table(het.file, head=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
het.outliers <- het[het$meanHet < mean(het$meanHet)-(3*sd(het$meanHet)) | het$meanHet > mean(het$meanHet)+(3*sd(het$meanHet)),]
het.output = paste(fileName,".het.outliers.ind",sep="")
write.table(het.outliers[,1:2],het.output,row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
if(!require(geneplotter)){
    source('http://bioconductor.org/biocLite.R')
    biocLite('geneplotter')
    }
colors  <- densCols(imiss$logF_MISS,het$meanHet)
imiss.het.pdf <- paste(fileName,".imiss.vs.het.pdf",sep="")
pdf(imiss.het.pdf)
    plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5),pch=20, xlab='Proportion of missing genotypes', ylab='Heterozygosity rate',axes=FALSE)
    axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=TRUE)
    axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
    abline(h=mean(het$meanHet)-(3*sd(het$meanHet)),col='RED',lty=2)
    abline(h=mean(het$meanHet)+(3*sd(het$meanHet)),col='RED',lty=2)
    abline(v=-1.522879, col='RED', lty=2)
dev.off()

# PCA plot and identify population outliers (>6*SD)
pca.file <- paste(fileName,".hapmap3r2.pca.evec",sep="")
data=read.table(pca.file,head=F,skip=1)
data$color <- 'black'
data$color[data$V12 %in% c('5','4')] <- 'PURPLE'
data$color[data$V12 %in% '3'] <- 'RED'
data$color[data$V12 %in% '6'] <- 'GREEN'
index=which(data$V12 %in% c('Case','Control'))
pca.pdf <- paste(fileName,".pca.pdf",sep="")
pdf(pca.pdf)
    plot(data$V2,data$V3,pch=20,xlab='PC1', ylab='PC2',col=data$color)
    points(data$V2[index],data$V3[index],pch=20,col='#00001110')
    legend('topright',c('CEU','CHB+JPT','YRI',fileName),col=c('RED','PURPLE','GREEN','BLACK'),pch=20)
    plot(data[,2:4],,col=data$color)
dev.off()
pca <- data[data$color %in% 'black',]
pc1.outliers <- pca[pca$V2 < (mean(data$V2[data$color %in% 'RED']) - 6*sd(data$V2[data$color %in% 'RED'])) | pca$V2 > (mean(data$V2[data$color %in% 'RED']) + 6*sd(data$V2[data$color %in% 'RED'])),]
pc2.outliers <- pca[pca$V3 < (mean(data$V3[data$color %in% 'RED']) - 6*sd(data$V3[data$color %in% 'RED'])) | pca$V3 > (mean(data$V3[data$color %in% 'RED']) + 6*sd(data$V3[data$color %in% 'RED'])),]
pca.outliers <- unique(rbind(pc1.outliers, pc2.outliers))
pca.outliers$ID <- gsub(':','\t',pca.outliers$V1)
pca.output <- paste(fileName, ".pca.outliers.ind", sep="")
write.table(pca.outliers$ID, pca.output, row.name=FALSE, col.names=F, quote=F, sep='\t')

# Extract the subject with higher missing rate for each related pair
related.file <- paste(fileName, ".related.ids", sep="")
ids <- read.table(related.file, head=TRUE, as.is=T)
imiss <- imiss[,c('FID','IID','F_MISS')]
names(imiss) <- c('FID2','IID2','F_MISS2')
ids <- merge(ids, imiss, sort=FALSE)
names(imiss) <- c('FID1','IID1','F_MISS1')
ids <- merge(ids,imiss,sort=F)
ids$FID <- ifelse(ids$F_MISS1 > ids$F_MISS2, ids$FID1, ids$FID2)
ids$IID <- ifelse(ids$F_MISS1 > ids$F_MISS2, ids$IID1, ids$IID2)
ids <- ids[,c('FID','IID','FID1','IID1','F_MISS1','FID2','IID2','F_MISS2','RT','EZ','Z0','Z1','Z2','PI_HAT')]
related.output <- paste(fileName, ".related.ind.info", sep="")
write.table(ids, related.output, row.names=F, quote=F, sep='\t')

