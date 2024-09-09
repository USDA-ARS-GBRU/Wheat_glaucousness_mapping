#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#rMVP GWAS
#https://github.com/xiaolei-lab/rMVP
#devtools:: install_github("xiaolei-lab/rMVP")
library(rMVP)
library(CMplot)
setwd("F:/Manuscript_wax_Daniela/2024_GWAS_Tshooting")
rm(list = ls())
pacman::p_load(snow,doSNOW,parallel,rrBLUP,dplyr,
               gaston,magrittr,lme4,lmerTest,emmeans,
               memisc,psych,sjPlot,tidyverse,tibble,
               popkin,factoextra,FactoMineR,popkin,
               BEDMatrix,vcfR,rMVP,CMplot,remotes,
               devtools,bwardr,gaston,dplyr,ggplot2,tidyr,tibble)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#_Step_1.1: Formatting/Intersect Pheno and geno data to make balance
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:->load raw BLUEs.csv
pheno1<-read.csv("WaxYearsBLUES.csv",sep=",",
                 head=TRUE,na.strings="NA")
pheno<-pheno1%>%distinct(Taxa,.keep_all=TRUE)#double check/clean there is not duplicated TAXA
head(pheno)
tail(pheno)
dim(pheno) #176

#:->load raw VCF file
gbs_mat<-read.vcf("2020_GAWN_SUNWHEAT_MDXN_6ST_postimp_filt.vcf.gz",
                  convert.chr=F)
dim(gbs_mat) #181 lines + 30519 SNP
class(gbs_mat)

#:->>Intersect Geno & Pheno data + write geno data as "genotypic_dat"
#to be uploaded & formatted below
genotypic_dat<-intersect(gbs_mat@ped$id,pheno$Taxa)#to filter VCF using the cvs file
genotypic_dat<-select.inds(gbs_mat, id %in% genotypic_dat)#attaching all VCF info
dim(genotypic_dat) #170 lines + 30519 SNP
class(genotypic_dat)
write.bed.matrix(genotypic_dat,"genotypic_dat",rds=NULL)
#NOTE: this will write 3 different files
# .bed , .bim and .fam, this will be used in Step_1.2.

#:->>Intersect Pheno and Geno data
#:->>>Extract gmat row names (VCF) to filter Pheno data sample IDs
GBS_as.matrix<-as.matrix(genotypic_dat) #make VFC a matrix to extract row_names
gbs_mat1<-GBS_as.matrix%>%
  as.data.frame(.)%>%
  tibble::rownames_to_column("Taxa")%>%
  dplyr::select(Taxa)
class(gbs_mat1)
head(gbs_mat1)
dim(gbs_mat1) #170
#write.csv(gbs_mat1,"./rMVP/gbs_mat1.csv")
#:->>>Filtering Pheno using Geno row names + load data as "phenotype_dat"
phenotype_dat<-inner_join(gbs_mat1,pheno,by=c("Taxa"))# ALL good
str(phenotype_dat)
dim(phenotype_dat) #170 + 3 variables
head(phenotype_dat,n=3)
tail(phenotype_dat,n=3)
#write.csv(SNB_all,"./rMVP/SNB_all.csv")

################################################################################
#_Step_1.2: Formatting Geno data
################################################################################
#rm(list=ls()) #clear R environment
#:->Full-featured function (Recommended)
MVP.Data(fileBed="genotypic_dat",
         filePhe=NULL,
         fileKin=TRUE,
         filePC=TRUE,       
         #priority="speed",
         #maxLine=10000,
         out="mvp.plink"
)
################################################################################
#_Step_2: Import formatted data (MVP. Data) from working directory
################################################################################
#phenotype_dat<-read.table("./rMVP/SNBblues.csv",head=TRUE) #already uploaded
genotypic_data<-attach.big.matrix("mvp.plink.geno.desc")
map_info<-read.table("mvp.plink.geno.map",head=TRUE)
#popstr<-read.table("mdp_population structure.txt",header=T,skip=1)
#Kinship<-MVP.K.VanRaden(genotypic_data,verbose = T)

####################################################################################################################################
#:
#:::::::::::::::::::::::::::::_Step_2: GWAS run:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#_2.1:->LOOP for several traits USING PCA as covariate ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for (i in 2:ncol(phenotype_dat)){
  imMVP<-MVP(
    phe=phenotype_dat[, c(1, i)],
    geno=genotypic_data,
    map=map_info,
    #K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=4,
    nPC.MLM=4,
    nPC.FarmCPU=4,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",#methods of variance components analysis, three methods are avaiblable, "BRENT", "EMMA", and "HE"
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM","MLM","FarmCPU")
  )
  gc()
}