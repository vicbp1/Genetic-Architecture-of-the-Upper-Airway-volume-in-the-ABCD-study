/local/chib/toconnor_grp/victor/softwares/julia/julia

using Interpolations, TextParse, DataFrames, MendelPlots, Missings, TrajGWAS, CSV, CategoricalArrays, LinearAlgebra
using Ipopt, WiSER
using Dates, DelimitedFiles
using Statistics
using Serialization

datadir="/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed_20250626/"

df = CSV.read("/local/chib/toconnor_grp/victor/IHC_projects/ABCD/Imputation_TOPMed/LongGWAS/TrajGWAS/20250525_ThreeTimesPoints_2143uninds_Sites_PCA.txt", DataFrame)

df.Sex = categorical(df.Sex)
df[!, :Age] = Float64.(df[!, :Age])
df[!, :Volume] = Float64.(df[!, :Volume])
df.Ethnicity = categorical(df.Ethnicity)
df.Site = categorical(df.Site)

for chr in 2
        trajgwas(@formula(Volume ~ 1 + Age + Sex + BMI_percentile + PC1 + PC2 + PC3 + PC4 + PC5 +PC6 +PC7 + PC8 +PC9 + PC10),
                @formula(Volume ~ 1),
                @formula(Volume ~ 1 + Age + Sex + BMI_percentile),
                :ids,
                df,
                datadir * "chr$chr/20250626_r2_0.8_3TP_2141_unrelated_chr$chr.assoc",
                pvalfile = datadir * "LongGWAS/TrajGWAS/20250630_Age_Sex_BMI_pcs1-10_SPA.chr$chr.pval.txt.gz",
                nullfile = datadir * "LongGWAS/TrajGWAS/20250630_Age_Sex_BMI_pcs1-10_SPA.chr$chr.null.txt",
                geneticformat = "VCF",
                vcftype = :DS, runs=10,usespa=true)
end
