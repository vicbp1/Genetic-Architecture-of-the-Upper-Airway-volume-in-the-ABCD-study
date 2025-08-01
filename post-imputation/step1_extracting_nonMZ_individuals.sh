#!/bin/bash
set -euo pipefail

chr=$1
path=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed_20250626
bcftools=/local/chib/toconnor_grp/victor/softwares/bcftools-1.19/bin/bcftools
plink2=/local/chib/toconnor_grp/victor/softwares/plink2
imputedvcf=$path/chr${chr}/chr${chr}.dose.vcf.gz
tempvcf0=$path/chr${chr}/20250626_r2_0.8_norm_all_chr${chr}.vcf.gz
tempvcf1=$path/chr${chr}/20250626_r2_0.8_updated_10204nonMZ_chr${chr}.vcf.gz
tempvcf2=$path/chr${chr}/20250626_r2_0.8_updated_2141unrelated_chr${chr}.vcf.gz
baseline_assoc=$path/chr${chr}/20250626_r2_0.8_Baseline_10204nonMZ_chr${chr}.assoc.vcf.gz
baseline_hsq=$path/chr${chr}/20250626_r2_0.8_Baseline_8802Unrelated_chr${chr}.hsq
trajgwasinput=$path/chr${chr}/20250626_r2_0.8_3TP_2141_unrelated_chr${chr}.assoc.vcf.gz
list=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/freeze_data/nonMZ_list_Baseline.txt 
related=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Population_Structure/20250503_Baseline_relatedness_0.0625_toRemove.txt
unrelated3TP=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/freeze_data/unrelated_list_3TP.ids
threads=4

## TOPMed Imputation results have multiallelic variants splitted as biallelic

$bcftools view -i 'R2>.8' $imputedvcf -Ou --threads ${threads} | $bcftools norm -m +any -Oz -o ${tempvcf0} --threads ${threads} ; $bcftools index ${tempvcf0}

## For Baseline: Input for Genesis and GCTA

$bcftools view -S $list ${tempvcf0} -Ou --threads ${threads} | $bcftools view -m2 -M2 -Ou --threads ${threads} \
| $bcftools +fill-tags -Oz -o ${tempvcf1} -- -t AC,MAF ; $bcftools index ${tempvcf1}

$bcftools view -i 'AC>20' ${tempvcf1} --threads ${threads} -Oz -o ${baseline_assoc} ; $bcftools index ${baseline_assoc}

$plink2 --vcf ${tempvcf1} --remove $related --maf 0.01 --make-bed --set-all-var-ids @:# --out ${baseline_hsq} --memory 39500


## For Three Points: Input for TrajGWAS

$bcftools view -S $unrelated3TP ${tempvcf0} -Ou --threads ${threads} | $bcftools view -m2 -M2 -Ou --threads ${threads} \
| $bcftools +fill-tags -Oz -o ${tempvcf2} -- -t AC,MAF ; $bcftools index ${tempvcf2}

$bcftools view -i 'AC>20' ${tempvcf2} --threads ${threads} -Oz -o ${trajgwasinput} ; $bcftools index ${trajgwasinput}


#  Cleanup (only if everything above succeeded) ===
rm -f ${tempvcf0} ${tempvcf0}.csi "$tempvcf1" "$tempvcf1.csi" "$tempvcf2" "$tempvcf2.csi"

