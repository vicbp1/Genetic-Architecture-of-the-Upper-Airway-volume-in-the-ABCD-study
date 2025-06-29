
### Running KING to identify MZ and to removed from PCA and association analyses

library("GENESIS")
library(gdsfmt)
library(SNPRelate)
library(SeqVarTools)
library(GWASTools)

## Running PCA
########################

path.samples<-"/Users/vborda/Documents/IHC/ABCD/plink_files/"
#name.file<-"20250515_Phenotypes_ThreeTimesPoints_GenotypingData"
name.file<-"20250515_Phenotypes_Baseline_GenotypingData"

bed.fn <- paste0(path.samples,name.file,".bed")
bim.fn <- paste0(path.samples,name.file,".bim")
fam.fn <- paste0(path.samples,name.file,".fam")

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, paste0(path.samples,name.file,".gds"))
snpgdsSummary(paste0(path.samples,name.file,".gds"))

(gdsfile <- snpgdsOpen(paste0(path.samples,name.file,".gds")))

ibd.robust  <- snpgdsIBDKING(gdsfile, num.thread=8)

path<-"/Users/vborda/Documents/IHC/ABCD/Population_Structure/"
save(ibd.robust,file=paste0(path,"20250515_Baseline_step2_KING_Results.RData"))
#save(ibd.robust,file=paste0(path,"20250515_ThreeTimesPoints_step2_KING_Results.RData"))
#load(file=paste0(path,"20250515_Baseline_step2_KING_Results.RData"))

KINGmat_11K<-kingToMatrix(ibd.robust)

## Getting MZ inds IDs

indices <- which(KINGmat_11K > 0.4 & row(KINGmat_11K) != col(KINGmat_11K), arr.ind = TRUE)

MZlist <- data.frame(
  Row_ID = rownames(KINGmat_11K)[indices[, 1]],  # Row names
  Column_ID = colnames(KINGmat_11K)[indices[, 2]],  # Column names
  Value = KINGmat_11K[indices]
)

fullids<-colnames(KINGmat_11K)

list_nonMZ <- fullids[!(fullids %in% MZlist$Row_ID)]

### LD Prunning

snpset_nonMZ <- snpgdsLDpruning(gdsfile, sample.id=list_nonMZ, method="corr", slide.max.bp=10e6, 
                                    ld.threshold=0.4, verbose=FALSE)
pruned_nonMZ <- unlist(snpset_nonMZ, use.names=FALSE)

mypcair_nonMZ <- pcair(gdsfile, sample.include=list_nonMZ,kinobj = KINGmat_nonMZ, divobj = KINGmat_nonMZ, kin.thresh=2^(-9/2), div.thresh=-2^(-9/2),
                           snp.include = pruned_nonMZ, maf=0.05,num.cores=8,eigen.cnt=50)

path<-"/Users/vborda/Documents/IHC/ABCD/Population_Structure/"
save(mypcair_nonMZ,file=paste0(path,"20250426_11k_step3_PCAir_nonMZ.RData"))

#### Running PCRelate

snpgdsClose(gdsfile) ### Closing since there is an issue with the gdsfile reading file

Genotypes <- GdsGenotypeReader(filename = paste0(path.samples,name.file,".gds"))
GenoData <- GenotypeData(Genotypes)

genoData <- GenotypeBlockIterator(GenoData, snpInclude=pruned_nonMZ)

mypcrelate_Baseline_nonMZ <- pcrelate(genoData, pcs = mypcair_nonMZ$vectors[,1:4], sample.include=list_nonMZ,
                       training.set = mypcair_nonMZ,
                       BPPARAM = BiocParallel::SerialParam())

save(mypcrelate_Baseline_nonMZ,file=paste0(path,"20250530_step3_PCRelate_10204ins.RData"))

### IBD probabilities

## Input for NATORA
IBDtable_baseline<-as.data.frame(mypcrelate_Baseline_nonMZ$kinBtwn)
write.table(IBDtable_baseline[,c(1,2,3)],file=paste0(path,"Population_Structure/20250530_10204inds_kinship.txt"),sep = '\t',row.names=FALSE,quote=FALSE,col.names = FALSE)
