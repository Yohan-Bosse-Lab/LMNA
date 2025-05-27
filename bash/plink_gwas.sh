#cd
cd data

#calculate Principal Components based on genotypic data
/home/renseb01/plink2 \
  --vcf LMNA_52_phased/chr_52.vcf \
  --pca 4 \
  --out LMNA_52_phased/pca_loadings

#Sort the output from above based on the alphabetic sample IDs.
(head -n 1 LMNA_52_phased/pca_loadings.eigenvec && tail -n +2 LMNA_52_phased/pca_loadings.eigenvec | sort -k1,1) >LMNA_52_phased/covariate.tsv

#add it to the covariates.
paste covariate.tsv severe.tsv >.LMNA_52_phased/phenotype.txt

#Do the glm in plink (no covariates)
#/home/renseb01/plink2 \
#  --vcf ../LMNA_52_phased/chr_52.vcf \
#  --pheno ../LMNA_52_phased/severe.tsv \
#  --pheno-name severity_recoded \
#  --glm allow-no-covars \
#  --out ../LMNA_52_phased/gwas_results

#glm with covariates
/home/renseb01/plink2 \
  --vcf ../LMNA_52_phased/chr_52.vcf \
  --pheno ../LMNA_52_phased/phenotype.txt \
  --pheno-name severity_recoded \
  --covar-name sex_recoded PC1 PC2 PC3 \
  --glm \
  --out ../LMNA_52_phased/gwas_results

