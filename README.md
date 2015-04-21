GWAS_QC.sh: a script that can be used to QC your GWAS data.
Please see this paper for more details: http://www.nature.com/nprot/journal/v5/n9/full/nprot.2010.116.html

```
Usage: GWAS_QC.sh <root name of the GWAS files> genome_build

If you have a gwas data set in plink binary format named as my_gwas.bed my_gwas.bim my_gwas.fam,
and the base pair positons in the bim file are based on the genome build hg18.
You can just use the following command to QC your data:

bash GWAS_QC.sh my_gwas hg18
```
