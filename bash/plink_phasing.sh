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

#filter MAF <0.05
vcftools --vcf YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples.vcf --maf 0.05 --out YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples_MAF005 --recode

#remove duplicates and filter unusable regions.
awk '!($1 == 1 && $2 < 2160299)' YB-P01_LORD_only_Plus_Steinberg_LMNA_148_samples_MAF005.recode.vcf >temp.vcf
awk '$1 != "0" && $1 != "X" && $1 != "XY" && $1 != "Y" && $1 != "MT" && $4 != "D" && $4 != "I" && $4 != "N"' temp.vcf >temp.vcf2
awk '!seen[$1, $2]++' temp.vcf2 | grep -iv 'ilmnseq' |  grep -iv 'Dup' >filtered.vcf

#phasing
beagle gt=filtered.vcf out=chr_phased

#Keep only chr1
gzip -d chr_phased.vcf
head -9 chr_phased.vcf >chr1_phased.vcf
awk '$1 == "1"' chr_phased.vcf >>chr1_phased.vcf
gzip chr1_phased.vcf

#cleanup
rm temp*


#THE REST is DONE via an Rscript locally
