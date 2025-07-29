## Estimating Heritability from Twin data

#source('https://vipbg.vcu.edu/vipbg/OpenMx2/software/getOpenMx.R')
library("OpenMx")
require(psych)

path.demo<-"/Users/vborda/Documents/IHC/ABCD/"
load(file=paste0(path.demo,"KinshipCoefficient_20250724_Baseline_MZtwins.RData"))
load(file=paste0(path.demo,"KinshipCoefficient_20250724_Baseline_DZtwins.RData"))
path.baseline.data<-"/Users/vborda/Documents/IHC/ABCD/"
load(file=paste0(path.baseline.data,"/Demographics/20250426_Phenotypes_Baseline.RData"))
MZtwins<-MZlist
rm(MZlist)
DZtwins$zygocity<-"DZ"
MZtwins$zygocity<-"MZ"

Twins<-rbind(DZtwins[,c("Column_ID","Row_ID","zygocity")],MZtwins[,c("Column_ID","Row_ID","zygocity")])

Twins.covariates<-merge(Twins,fulldata_baseline[,c("ids2","Volume","Age","BMI_percentile","Sex","Ethnicity")],by.x = "Column_ID",by.y="ids2")
colnames(Twins.covariates)<-c("ids1","ids2","zygocity","Volume1","Age1","BMI_percentile1","Sex1","Ethnicity1")
Twins.covariates<-merge(Twins.covariates,fulldata_baseline[,c("ids2","Volume","Age","BMI_percentile","Sex","Ethnicity")],by.x = "ids2",by.y="ids2")
colnames(Twins.covariates)[9:13]<-c("Volume2","Age2","BMI_percentile2","Sex2","Ethnicity2")


Twins.covariates$key <- apply(Twins.covariates[, c("ids1", "ids2")], 1, function(x) paste(sort(x), collapse = "_"))
Twins.covariates_unique <- Twins.covariates[!duplicated(Twins.covariates$key), ]
Twins.covariates_unique$key <- NULL
Twins.covariates<-Twins.covariates_unique

Twins.covariates$Volume1_z <- scale(Twins.covariates$Volume1)
Twins.covariates$Volume2_z <- scale(Twins.covariates$Volume2)

#############################################################
####### Heritability calculation adjusting for covariates    
#############################################################
# MZ and DZ data
nv<-1
nt<-2
ntv<-nv*nt



Twins.covariates$Volume1_z_reg <- resid(lm(data = Twins.covariates, Volume1_z ~ Age1 + Sex1 + BMI_percentile1 + Ethnicity1, na.action = na.exclude))
Twins.covariates$Volume2_z_reg <- resid(lm(data = Twins.covariates, Volume2_z ~ Age2 + Sex2 + BMI_percentile2 + Ethnicity2, na.action = na.exclude))

selVars <- c("Volume1_z_reg", "Volume2_z_reg")
mzData <- subset(Twins.covariates, zygocity == "MZ", select = selVars)
dzData <- subset(Twins.covariates, zygocity == "DZ", select = selVars)

# For St_a c and e (additive genetic)
St_a <- St_c <- St_e <- sd(Twins.covariates$Volume_resid, na.rm=T) * 0.3

### Fit ACE Model with RawData and Matrices Input
# ------------------------------------------------------
pathA <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=St_a,
                   labels="a11", name="a")
#pathC <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=St_c,
#                   labels="c11", name="c")
# Set C to zero and don't estimate it
pathC <- mxMatrix("Lower", nv, nv, free=FALSE, values=0, labels="c11", name="c")

pathE <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=St_e,
                   labels="e11", name="e")
covA <- mxAlgebra( expression=a %*% t(a), name="A")
covC <- mxAlgebra( expression=c %*% t(c), name="C")
covE <- mxAlgebra(expression=e %*% t(e), name="E")
covP <- mxAlgebra( expression=A+C+E, name="V")
matI <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

covMZ <- mxAlgebra(rbind(cbind(A+C+E, A+C),
                         cbind(A+C, A+C+E)), name="expCovMZ")

stpatha <- mxAlgebra( iSD %*% a, name="sta")
stpathc <- mxAlgebra( iSD %*% c, name="stc")
stpathe <- mxAlgebra( iSD %*% e, name="ste")
stcovA <- mxAlgebra( A/V, name="stA")
stcovC <- mxAlgebra( C/V, name="stC")
stcovE <- mxAlgebra( E/V, name="stE")


Stmean <- vech(mean(mzData$Volume1_z_reg,na.rm=T))

meanG <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=Stmean, label=c
                   ("mean","mean"), name="expMean" )
covMZ <- mxAlgebra( expression= rbind(cbind(A+C+E,A+C),
                                      cbind(A+C,A+C+E)), name="expCovMZ")
covDZ <- mxAlgebra( expression= rbind(cbind(A+C+E , 0.5%x%A+C),
                                      cbind(0.5%x%A+C,A+C+E)), name="expCovDZ")

dataMZ <- mxData( observed=mzData, type="raw")
dataDZ <- mxData( observed=dzData, type="raw")
expMZ <- mxExpectationNormal(covariance="expCovMZ",
                             means="expMean", dimnames=selVars)

expDZ <- mxExpectationNormal( covariance="expCovDZ",
                              means="expMean", dimnames=selVars )
fitfun <- mxFitFunctionML()

pars <- list(pathA, pathC, pathE, covA, covC, covE, covP, matI, invSD, stcovA,
             stcovC, stcovE, stpatha, stpathc, stpathe)
modelMZ <- mxModel( pars, meanG, covMZ, dataMZ, expMZ, fitfun, name="MZ")
modelDZ <- mxModel( pars, meanG, covDZ, dataDZ, expDZ, fitfun, name="DZ")
multi <- mxFitFunctionMultigroup( c("MZ","DZ"))
CIs <- mxCI( c('stA','stC','stE'))
modelAce <- mxModel( "ACE", pars, modelMZ, modelDZ, multi, CIs)
fitAce <- mxTryHard(modelAce, intervals=TRUE)

summary(fitAce)
