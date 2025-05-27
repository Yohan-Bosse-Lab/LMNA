#!/bin/bash

#original data is here:
#mkdir data/LMNA_148_samples
#cp /home/renseb01/Documents/mendel_Christian_Steinberg/YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples/PLINK_200225_0230/* data/LMNA_148_samples/.

#go to data/
cd data/LMNA_148_samples

#sort and  make pgen
/home/renseb01/plink2 --ped YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples.ped --map YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples.map --make-pgen --out YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples --allow-extra-chr --sort-vars

#make the fam file
/home/renseb01/plink2 --pfile YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples --make-bed --out YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples --allow-extra-chr

#make a .vcf file
/home/renseb01/plink2 --bfile YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples --recode vcf --out YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples

#remove duplicates, keep only chromosome1, and filter unusable regions.
#filter MAF <0.05
vcftools --vcf YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples.vcf --maf 0.05 --out YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples_MAF005 --recode

head -33 YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples_MAF005.recode.vcf >header

awk '!($1 == 1 && $2 < 2160299)' YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples_MAF005.recode.vcf >temp.vcf
awk '$1 != "X"' temp.vcf >temp.vcf2
awk '$1 != "XY"' temp.vcf2 >temp.vcf3
awk '$1 != "Y"' temp.vcf3 >temp.vcf4
awk '$1 != "MT"' temp.vcf4 >temp.vcf5
awk '!seen[$1, $2]++' temp.vcf5 >temp.vcf6
awk '$4 != "D"' temp.vcf6 >temp.vcf7
awk '$4 != "I"' temp.vcf7 >temp.vcf8
awk '$4 != "N"' temp.vcf8 >temp.vcf9
grep -iv 'ilmnseq' temp.vcf9 >temp.vcf10
grep -iv 'Dup' temp.vcf10 >temp.vcf11
cat temp.vcf11 >>header
mv header filtered.vcf

#phasing
beagle gt=filtered.vcf out=chr_phased

#reformat for chr1
zcat chr_phased.vcf.gz | head -9 >chr1_phased.vcf
zcat chr_phased.vcf.gz | awk '$1 == "1"' >>chr1_phased.vcf
gzip chr1_phased.vcf

#cleanup
rm temp*
rm header

#THE REST is DONE via an Rscript locally
