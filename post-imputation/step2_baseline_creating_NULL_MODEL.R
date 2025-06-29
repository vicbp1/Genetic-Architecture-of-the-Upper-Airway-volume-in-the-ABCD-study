library("GENESIS")
library(gdsfmt)
library(SNPRelate)
library(SeqVarTools)
library(GWASTools)

### Extracting Phenotypes for Baseline individuals 

path<-"/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Demographic_Data/"
path2<-"/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Population_Structure/"

load(file=paste0(path,"20250530_nonMZ_Baseline_10204inds_covar_pcs_final.RData"))
load(file=paste0(path,"20250530_MatrixSites_Baseline.RData"))
load(file=paste0(path2,"20250530_step3_PCRelate_10204ins.RData"))

myGRM <- pcrelateToMatrix(mypcrelate_Baseline_nonMZ)
myGRMrun<-as.matrix(myGRM)
scanAnnot <- ScanAnnotationDataFrame(nonMZ_baseline_genotyping_equalSex_PCA)

# Creating null model with two random variables

nullmod_20250530_10204inds_ethnicity_2RV <- fitNullModel(scanAnnot,outcome="Volume",covars=c("Age", "Sex" ,"BMI_percentile",
                                                                                   "Ethnicity", "PC1" , "PC2" , "PC3" , "PC4",
                                                                                   "PC5", "PC6", "PC7","PC8", "PC9", "PC10"),
                                               cov.mat = list(nonMZ_baseline_site_matrix,myGRMrun), family = "gaussian")

save(nullmod_20250530_10204inds_ethnicity_2RV,file="/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed/GWAS_Baseline/inputs/20250530_NULLMODEL_10204inds_ethnicity_2RV.RData")

#Testing Heritability

varCompCI(nullmod.height, prop = TRUE)
