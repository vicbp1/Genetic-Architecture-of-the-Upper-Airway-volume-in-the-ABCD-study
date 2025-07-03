#!/bin/bash
chr=$1
plink2=/local/chib/toconnor_grp/victor/softwares/plink2 
gcta=/local/chib/toconnor_grp/victor/softwares/gcta/gcta-1.94.1
path=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed_20250626/chr${chr}
unrelated=20250626_r2_0.8_Baseline_8802Unrelated_chr${chr}.hsq
outpath=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed_20250626/Heritability/chr${chr}
pathenv=/local/chib/toconnor_grp/victor/miniconda3/envs/RunningJasmine/envs/GWAS/bin
resources=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed/scripts
outfix=biallelic_unrelated_r2_0.8_chr${chr}
threads=5

$gcta --bfile $path/${unrelated} --ld-score-region 200 --out $outpath/unrelated_step1_r0.8.chr$chr --thread-num ${threads}

$pathenv/Rscript ${resources}/ldscore_bin.R $outpath/unrelated_step1_r0.8.chr$chr $chr $outpath/ ${outfix}

#estimate GRM stratified by LD score 
for group in $(seq 1 4); do
        $gcta --bfile $path/${unrelated} --extract $outpath/${outfix}.snp_group${group}.txt --make-grm --out $outpath/${outfix}.test_group${group} --thread-num ${threads}
done
