#!/local/chib/toconnor_grp/victor/miniconda3/envs/RunningJasmine/envs/GWAS/bin/Rscript
## To run
##
path=/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed
path2=$path/GWAS_Baseline
## for chr in 6; do Rscript 20250603_GWAS_Baseline.R $path2/inputs/20250530_NULLMODEL_10204inds_ethnicity_2RV.RData $path/chr${chr}/20250530_r2_0.3_Baseline_10204inds_chr${chr}.vcf.gz $path2/inputs/20250603_chr${chr}_baseline_10204inds.gds $path2/outputs/20250603_chr${chr}_10910inds_pcs_age_bmi_sex_ethnicity_2RV.txt ; done

args = commandArgs(trailingOnly=TRUE)

library("GENESIS")
library(SeqVarTools)
library(GWASTools)

## loading Models

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

NULLMODEL <- loadRData(args[1])

vcffile<-args[2]
gdsfileV<-paste0(args[3])
ncores=8

seqVCF2GDS(vcffile, gdsfileV, fmt.import="DS",
           storage.option="LZMA_RA", verbose=FALSE,
           parallel = ncores)
gds<-seqOpen(gdsfileV)

seqData <- SeqVarData(gds)

ABCD_genoData <- SeqVarBlockIterator(seqData, verbose=FALSE)

## Running SNP-Phenotype Association Test

assoc <- assocTestSingle(ABCD_genoData, null.model = NULLMODEL,
                         BPPARAM = BiocParallel::SerialParam(), imputed=TRUE)

write.table(assoc,args[4], sep='\t', col.names = TRUE, row.names = F, quote = F)
