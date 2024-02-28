# LA RIL pop - QTL Mapping - Raleigh 2021 - Wax/Glaucousness Phenotypes - Leaf + Heads

#options(max.print=999999)
#memory.limit(size = NA)
library("qtl")
library("magrittr")
library("ASMap")
#library("dplyr")

LAmap<-read.cross(format="csv",file="D:/wax/Field2021/LAGenMap_26LG+waxKASP_2021.03.01+RAL21+RAL20+fieldwax+avgs+K20.csv",estimate.map=FALSE,
                  na.strings=c(".","-","NA"),genotypes=c("A","B"),crosstype="riself")%>%jittermap(.)%>%calc.genoprob(map.function="kosambi")
summary(LAmap)


#### R20 LEAF BOTTOM ####
#column 6

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR20leafB<-scanone(LAmaprobR21,pheno.col=6,n.perm=5000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR20leafB, alpha=c(0.20,0.05))
#LOD thresholds (5000 permutations)
#lod
#20% 2.76
#5%  3.42

#save LOD result
setwd("D:/wax/Field2021/QTL-results")
save(out.permLAR20leafB, file = "R20-LA-wax-leafB-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R20-LA-wax-leafB-LODperm")


# SIM
SIMLAR20leafB<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR20leafB<-scanone(SIMLAR20leafB,pheno.col=6,model="normal",method="em")
summary(out.SIMLAR20leafB,alpha=0.05,perms=out.permLAR20leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#c2D.loc37      2D  37     29      57  5.46 8e-04
#S3A_603900770  3A 163    162     165 21.33 0e+00
#c4B.loc123     4B 123     90     136  4.38 6e-03

# CIM
#window size 6; 4 covariate markers
out.CIMLAR20leafB<-cim(SIMLAR20leafB,pheno.col=6, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR20leafB,alpha=0.05,perms=out.permLAR20leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod   pval
#c2B.loc11      2B  11.0   7.37    14.0  3.61 0.0334
#S2D_22444654   2D  47.2  45.00    51.8  5.70 0.0006
#S3A_603900770  3A 163.1 162.00   165.0 24.25 0.0000

#save QTL results
setwd("D:/wax/Field2021/QTL-results")
write.table(out.SIMLAR20leafB,file="LA_R20_waxleafB_SIM.csv", sep = ",")
write.table(out.CIMLAR20leafB,file="LA_R20_waxleafB_CIM.csv", sep = ",")


#### 4: Plots

#SIM plots
plot(out.SIMLAR20leafB, main = "Standard Interval Mapping - LA R20 wax leaf bottom")
add.threshold(out.SIMLAR20leafB,perms=out.permLAR20leafB,alpha=0.05, col="green")

#CIM plots
plot(out.CIMLAR20leafB, main = "Composite Interval Mapping - LA R20 wax leaf bottom")
add.threshold(out.CIMLAR20leafB,perms=out.permLAR20leafB,alpha=0.05, col="blue")

#significant QTL LGs only 
plot(out.SIMLAR20leafB,out.CIMLAR20leafB, chr=c("2B", "2D", "3A", "4B"), col=c("cyan", "blue"), main="LA pop Raleigh 2020 wax leaf bottom")
add.threshold(out.SIMLAR20leafB,perms=out.permLAR20leafB,alpha=0.05, col="orange")

#### R20 LEAF BOTTOM WAX - MODEL PVE AND EPISTASIS ####


#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#2D
R20.2D.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2D",pos=37,qtl.name="c2D.loc37",what="prob")
R20.2D.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=6,qtl=R20.2D.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.2D.SIM.spike.QTL.Model)

#3A
R20.3A.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163,qtl.name="S3A_603900770",what="prob")
R20.3A.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=6,qtl=R20.3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.3A.SIM.spike.QTL.Model)

#4B
R20.4B.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=123,qtl.name="c4B.loc123",what="prob")
R20.4B.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=6,qtl=R20.4B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.4B.SIM.spike.QTL.Model)


#CIM:

#2B
R20.2B.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=11,qtl.name="c2B.loc11",what="prob")
R20.2B.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=6,qtl=R20.2B.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.2B.CIM.spike.QTL.Model)

#2D
R20.2D.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2D",pos=47.2,qtl.name="S2D_22444654",what="prob")
R20.2D.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=6,qtl=R20.2D.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.2D.CIM.spike.QTL.Model)


## Fit a multiple-QTL model 
R20spikeQTLs<- makeqtl(LAmaprobR21,chr=c("2D","3A","4B"),pos=c(37,163,123),
                       qtl.name=c("R20.2D.SIM.spike.QTL","R20.3A.SIM.spike.QTL","R20.4B.SIM.spike.QTL"),
                       what="prob")
R20spikeQTL.model<- fitqtl(LAmaprobR21,pheno.col=2,qtl=R20spikeQTLs,formula=y~Q1*Q2*Q3,method="hk")
summary(R20spikeQTL.model)

#### No Significant epistasis 



#### R20 SPIKE WAX ####
#column 5

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR20spike<-scanone(LAmaprobR21,pheno.col=5,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR20spike, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.67
#5%  3.45

#save LOD result
setwd("D:/wax/Field2021/QTL-results/")
save(out.permLAR20spike, file = "R20-LA-wax-spike-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R20-LA-wax-spike-LODperm")


# SIM
SIMLAR20spike<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR20spike<-scanone(SIMLAR20spike,pheno.col=5,model="normal",method="em")
summary(out.SIMLAR20spike,alpha=0.05,perms=out.permLAR20spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#K1RS.6110     1B.1   0.0    0.0     3.0 15.64 0.000
#S2A_734831175   2A  69.7   67.4    75.0 17.58 0.000
#S3A_603900770   3A 163.1  156.0   170.4  3.71 0.023
#c4D.loc74       4D  74.0   47.0    80.8  5.24 0.002

# CIM
#window size 6; 4 covariate markers
out.CIMLAR20spike<-cim(SIMLAR20spike,pheno.col=5, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR20spike,alpha=0.05,perms=out.permLAR20spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod   pval
#K1RS.6110     1B.1   0.0    0.0     3.0 20.36 0.000
#S2A_734831175   2A  69.7   67.4    71.4 23.79 0.000
#c3A.loc163      3A 163.0  160.0   167.0  7.90 0.000
#c4D.loc62       4D  62.0   61.0    68.0  5.80 0.000
#c6B.loc183      6B 183.0  178.0   187.7  4.31 0.004
#c7A.2.loc80   7A.2  80.0   73.0    85.7  3.48 0.045

#save QTL results
setwd("D:/wax/Field2021/QTL-results")
write.table(out.SIMLAR20spike,file="LA_R20_waxspike_SIM.csv", sep = ",")
write.table(out.CIMLAR20spike,file="LA_R20_waxspike_CIM.csv", sep = ",")




#### R20 SPIKE WAX - MODEL PVE AND EPISTASIS ####


#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#1B.1
R20.1B.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="1B.1",pos=0,qtl.name="K1RS.6110",what="prob")
R20.1B.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.1B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.1B.SIM.spike.QTL.Model)

#2A
R20.2A.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=69.7,qtl.name="S2A_734831175",what="prob")
R20.2A.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.2A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.2A.SIM.spike.QTL.Model)

#3A
R20.3A.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163.1,qtl.name="S3A_603900770",what="prob")
R20.3A.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.3A.SIM.spike.QTL.Model)

#4D
R20.4D.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=74,qtl.name="c4D.loc74",what="prob")
R20.4D.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.4D.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.4D.SIM.spike.QTL.Model)



# CIM:

#3A
R20.3A.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163,qtl.name="c3A.loc163",what="prob")
R20.3A.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.3A.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.3A.CIM.spike.QTL.Model)

#4D
R20.4D.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=62,qtl.name="c4D.loc62",what="prob")
R20.4D.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.4D.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.4D.CIM.spike.QTL.Model)

#6B
R20.6B.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=183,qtl.name="c6B.loc183",what="prob")
R20.6B.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=5,qtl=R20.6B.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20.6B.CIM.spike.QTL.Model)


## Fit a multiple-QTL model 
R20spikeQTLs<- makeqtl(LAmaprobR21,chr=c("1B.1","2A","3A","4D"),pos=c(0,69.7,163.1,74),
                       qtl.name=c("R20.1B.SIM.spike.QTL","R20.2A.SIM.spike.QTL","R20.3A.SIM.spike.QTL","R20.4D.SIM.spike.QTL"),
                       what="prob")
R20spikeQTL.model<- fitqtl(LAmaprobR21,pheno.col=5,qtl=R20spikeQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(R20spikeQTL.model)





#### R21 LEAF BOTTOM ####
#column 3

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR21leafB<-scanone(LAmaprobR21,pheno.col=3,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR21leafB, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.78
#5%  3.42

#save LOD result
setwd("D:/wax/Field2021/QTL-results")
save(out.permLAR21leafB, file = "R21-LA-wax-leafB-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R21-LA-wax-leafB-LODperm")


# SIM
SIMLAR21leafB<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR21leafB<-scanone(SIMLAR21leafB,pheno.col=3,model="normal",method="em")
summary(out.SIMLAR21leafB,alpha=0.05,perms=out.permLAR21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#S2A_737386570   2A  73.5   35.5     105  3.63 0.032
#c3A.loc163      3A 163.0  162.0     165 19.53 0.000
#S4B_596346774   4B 109.2   86.0     138  4.39 0.005
#c5B.2.loc19   5B.2  19.0    6.0      34  4.02 0.012
#S6B_561009143   6B 133.7  114.0     145  3.43 0.050

# CIM
#window size 6; 4 covariate markers
out.CIMLAR21leafB<-cim(SIMLAR21leafB,pheno.col=3, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR21leafB,alpha=0.05,perms=out.permLAR21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod   pval
#S2A_737386570  2A  73.5     64      79  4.70 0.002
#c3A.loc163     3A 163.0    162     165 25.19 0.000
#S4B_595211656  4B 108.5    105     112  5.04 0.001
#S6B_558318258  6B 134.4    131     138  5.43 0.000

#save QTL results
setwd("D:/wax/Field2021/QTL-results")
write.table(out.SIMLAR21leafB,file="LA_R21_waxleafB_SIM.csv", sep = ",")
write.table(out.CIMLAR21leafB,file="LA_R21_waxleafB_CIM.csv", sep = ",")


#### R21 LEAF WAX BOTTOM - MODEL PVE AND EPISTASIS ####




#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#2A
R21.2A.SIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=73.5,qtl.name="S2A_737386570",what="prob")
R21.2A.SIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.2A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.2A.SIM.leaf.QTL.Model)

#3A
R21.3A.SIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163,qtl.name="c3A.loc163",what="prob")
R21.3A.SIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.3A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.SIM.leaf.QTL.Model)

#4B
R21.3A.SIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=109.2,qtl.name="S4B_596346774",what="prob")
R21.3A.SIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.3A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.SIM.leaf.QTL.Model)

#5B
R21.5B.SIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=19,qtl.name="c5B.2.loc19",what="prob")
R21.5B.SIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.5B.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.5B.SIM.leaf.QTL.Model)

#6B
R21.6B.SIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=133.7,qtl.name="S6B_561009143",what="prob")
R21.6B.SIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.6B.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.6B.SIM.leaf.QTL.Model)


# CIM:

#4B
R21.4B.CIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=108.5,qtl.name="S4B_595211656",what="prob")
R21.4B.CIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.4B.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4B.CIM.leaf.QTL.Model)

#6B
R21.6B.CIM.leaf.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=134.4,qtl.name="S6B_558318258",what="prob")
R21.6B.CIM.leaf.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=3,qtl=R21.6B.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.6B.CIM.leaf.QTL.Model)


## Fit a multiple-QTL model 
R21leafQTLs<- makeqtl(LAmaprobR21,chr=c("2A","3A","4B","5B.2"),pos=c(73.5,163,109.2,19),
                       qtl.name=c("R21.2A.SIM.leaf.QTL","R21.3A.SIM.leaf.QTL","R21.4B.SIM.leaf.QTL","R21.5B.SIM.leaf.QTL"),
                       what="prob")
R21leafQTL.model<- fitqtl(LAmaprobR21,pheno.col=3,qtl=R21leafQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(R21leafQTL.model)

#### No Significant epistasis 



#### R21 SPIKE WAX  ####
#column 2

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR21spike<-scanone(LAmaprobR21,pheno.col=2,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR21spike, alpha=c(0.20,0.05))
#LOD thresholds (5000 permutations)
#lod
#20% 2.72
#5%  3.48

#save LOD result
setwd("D:/wax/Field2021/QTL-results")
save(out.permLAR21spike, file = "R21-LA-wax-spike-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R21-LA-wax-spike-LODperm")


# SIM
SIMLAR21spike<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR21spike<-scanone(SIMLAR20leafB,pheno.col=2,model="normal",method="em")
summary(out.SIMLAR21spike,alpha=0.05,perms=out.permLAR21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#c1A.loc77       1A 77.0   67.0    91.6  4.22 0.013
#c1B.1.loc3    1B.1  3.0    0.0     7.0 15.03 0.000
#c2A.loc74       2A 74.0   71.4    79.0 12.25 0.000
#S4B_551081152   4B 96.1   73.0   124.0  4.56 0.009
#c4D.loc77       4D 77.0   71.0    80.8 10.43 0.000

# CIM
#window size 6; 4 covariate markers
out.CIMLAR21spike<-cim(SIMLAR21spike,pheno.col=2, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR21spike,alpha=0.05,perms=out.permLAR21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod   pval
#S1A_33219502    1A  75.0     71    78.0  6.34 0.002
#c1B.1.loc1    1B.1   1.0      0     5.0 18.31 0.000
#S2A_737386570   2A  73.5     72    75.7 14.12 0.000
#S3A_603900770   3A 163.1    158   168.0  3.62 0.037
#c4B.loc121      4B 121.0     87   125.7  5.45 0.004
#c4D.loc80       4D  80.0     77    80.8 14.20 0.000
#c6B.loc186      6B 186.0    179   187.7  6.11 0.002

#save QTL results
setwd("D:/wax/Field2021/QTL-results")
write.table(out.SIMLAR21spike,file="LA_R21_waxspike_SIM.csv", sep = ",")
write.table(out.CIMLAR21spike,file="LA_R21_waxspike_CIM.csv", sep = ",")



#### R21 SPIKE WAX - MODEL PVE AND EPISTASIS ####


#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#1A
R21.1A.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="1A",pos=77,qtl.name="S1A_33236936",what="prob")
R21.1A.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.1A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1A.SIM.spike.QTL.Model)

#1B
R21.1B.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="1B.1",pos=3,qtl.name="c1B.1.loc3",what="prob")
R21.1B.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.1B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1B.SIM.spike.QTL.Model)

#2A
R21.2A.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=74,qtl.name="c2A.loc74",what="prob")
R21.2A.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.2A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.2A.SIM.spike.QTL.Model)

#4B
R21.4B.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=96.1,qtl.name="S4B_551081152",what="prob")
R21.4B.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.4B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4B.SIM.spike.QTL.Model)

#4D
R21.4D.SIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=77,qtl.name="c4D.loc77",what="prob")
R21.4D.SIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.4D.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4D.SIM.spike.QTL.Model)


#CIM:

#1A
R21.1A.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="1A",pos=75,qtl.name="S1A_33219502",what="prob")
R21.1A.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.1A.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1A.cIM.spike.QTL.Model)

#1B
R21.1B.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="1B.1",pos=1,qtl.name="c1B.1.loc1",what="prob")
R21.1B.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.1B.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1B.cIM.spike.QTL.Model)

#2A
R21.2A.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=73.5,qtl.name="S2A_737386570",what="prob")
R21.2A.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.2A.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.2A.cIM.spike.QTL.Model)


#3A
R21.3A.CIM.spike.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163.1,qtl.name="S3A_603900770",what="prob")
R21.3A.CIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.3A.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.CIM.spike.QTL.Model)

#4B
R21.4B.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=121,qtl.name="c4B.loc121",what="prob")
R21.4B.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.4B.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4B.cIM.spike.QTL.Model)

#4D
R21.4D.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=80,qtl.name="c4D.loc80",what="prob")
R21.4D.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.4D.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4D.cIM.spike.QTL.Model)

#6B
R21.6B.cIM.spike.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=186,qtl.name="c6B.loc186",what="prob")
R21.6B.cIM.spike.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=2,qtl=R21.6B.cIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.6B.cIM.spike.QTL.Model)


## Fit a multiple-QTL model 
R21spikeQTLs<- makeqtl(LAmaprobR21,chr=c("1A","1B.1","2A","4B","4D"),pos=c(77,3,74,96.1,77),
                       qtl.name=c("R21.1A.cIM.spike.QTL","R21.1B.cIM.spike.QTL","R21.2A.cIM.spike.QTL", "R21.4B.SIM.spike.QTL", "R21.4D.SIM.spike.QTL"),
                       what="prob")
R21spikeQTL.model<- fitqtl(LAmaprobR21,pheno.col=2,qtl=R21spikeQTLs,formula=y~Q1*Q2*Q3*Q4*Q5,method="hk")
summary(R21spikeQTL.model)

#### No Significant epistasis 










#### R21 + R20 LEAF BOTTOM average ####
#pheno column 9

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permleafBavg<-scanone(LAmaprobR21,pheno.col=9,n.perm=5000, verbose=TRUE,method="em")#,model="normal")
summary(out.permleafBavg, alpha=c(0.20,0.05))
#LOD thresholds (5000 permutations)
#lod
#20% 2.76
#5%  3.47

#save LOD result
setwd("D:/wax/Field2021/QTL-results")
save(out.permleafBavg, file = "R20-R21-LA-wax-leafB-avg-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-leafB-avg-LODperm")
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-leafB-avg-SIM-results")
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-leafB-avg-CIM-results")



# SIM
SIMLARavgleafB<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLARavgleafB<-scanone(SIMLARavgleafB,pheno.col=9,model="normal",method="em")
summary(out.SIMLARavgleafB,alpha=0.05,perms=out.permleafBavg,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
save(out.SIMLARavgleafB, file = "R20-R21-LA-wax-leafB-avg-SIM-results")
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod   pval
#c2D.loc38       2D  38   29.0      58  4.01 0.0120
#c3A.loc163      3A 163  162.0     165 26.51 0.0000
#S4B_635945715   4B 124   88.8     136  5.40 0.0002
#c5B.2.loc19   5B.2  19    0.0      34  3.93 0.01
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-leafB-avg-SIM-results")

# CIM
#window size 6; 4 covariate markers
out.CIMLARavgleafB<-cim(SIMLARavgleafB,pheno.col=9, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLARavgleafB,alpha=0.05,perms=out.permleafBavg,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
save(out.CIMLARavgleafB, file = "R20-R21-LA-wax-leafB-avg-CIM-results")
#    Length Class      Mode
#lod 6      data.frame list
#chr    pos ci.low ci.high   lod   pval
#S2A_737386570  2A  73.50     62  105.40  5.48 0.0002
#S2B_9713222    2B   1.31      0    4.54  5.47 0.0002
#c3A.loc163     3A 163.00    162  165.00 33.82 0.0000
#S4B_595211656  4B 108.45    105  112.00  4.92 0.0012
#S6B_558318258  6B 134.42    114  144.00  4.13 0.0080
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-leafB-avg-CIM-results")


#save QTL results
setwd("D:/wax/Field2021/QTL-results")
write.table(out.SIMLARavgleafB,file="LA_20R21_avg_waxleafB_SIM.csv", sep = ",")
write.table(out.CIMLARavgleafB,file="LA_20R21_avg_waxleafB_CIM.csv", sep = ",")


#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

LBM.2D.QTL<-makeqtl(LAmaprobR21,chr="2D",pos=38,qtl.name="c2D.loc38",what="prob")
LBM.2D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.2D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.2D.QTL.Model)

LBM.3A.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163,qtl.name="c3A.loc163",what="prob")
LBM.3A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.3A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.3A.QTL.Model)

LBM.4B.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=124,qtl.name="S4B_635945715",what="prob")
LBM.4B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.4B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.4B.QTL.Model)

LBM.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=19,qtl.name="c5B.2.loc19",what="prob")
LBM.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.5B.2.QTL.Model)

## Fit a multiple-QTL model 
LBMQTLs<- makeqtl(LAmaprobR21,chr=c("2D","3A", "4B","5B.2"),pos=c(38,163,124,19),
                    qtl.name=c("LBM.2D.QTL","LBM.3A.QTL", "LBM.4B.QTL", "LBM.5B.2.QTL"),
                    what="prob")
LBMMQTL.model<- fitqtl(LAmaprobR21,pheno.col=9,qtl=LBMQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(LBMMQTL.model)

# b) out.CIM
## Fit single QTL model
LBM.2A.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=73.5,qtl.name="S2A_737386570",what="prob")
LBM.2A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.2A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.2A.QTL.Model)

LBM.2B.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=1.31,qtl.name="S2B_9713222",what="prob")
LBM.2B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.2B.QTL.Model)

LBM.3A.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=163,qtl.name="c3A.loc163",what="prob")
LBM.3A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.3A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.3A.QTL.Model)

CIM.LBM.4B.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=108.45,qtl.name="S4B_595211656",what="prob")
CIM.LBM.4B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=CIM.LBM.4B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIM.LBM.4B.QTL.Model)

LBM.6B.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=134.42,qtl.name="S6B_558318258",what="prob")
LBM.6B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=9,qtl=LBM.6B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(LBM.6B.QTL.Model)


## Fit a multiple-QTL model 
LBMCIMMQTLs<- makeqtl(LAmaprobR21,chr=c("2A", "2B", "3A", "4B", "6B"),pos=c(73.5,1.31,163,108.45,134.42),
                       qtl.name=c("LBM.2A.QTL","LBM.2B.QTL", "LBM.3A.QTL", "CIM.LBM.4B.QTL", "LBM.6B.QTL"),
                       what="prob")
LBMCIMMQTL.model<- fitqtl(LAmaprobR21,pheno.col=9,qtl=LBMCIMMQTLs,formula=y~Q1*Q2*Q3*Q4*Q5,method="hk")
summary(LBMCIMMQTL.model)



#### R21 + R20 SPIKE average ####
#pheno column 8

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permspikeavg<-scanone(LAmaprobR21,pheno.col=8,n.perm=5000, verbose=TRUE,method="em")#,model="normal")
summary(out.permspikeavg, alpha=c(0.20,0.05))
#LOD thresholds (5000 permutations)
#lod
#20% 2.75
#5%  3.43

#save LOD result
save(out.permspikeavg, file = "D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-LODperm")
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-SIM-results")
load("D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-CIM-results")

# SIM
SIMLAavgspike<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAavgspike<-scanone(SIMLAavgspike,pheno.col=8,model="normal",method="em")
summary(out.SIMLAavgspike,alpha=0.05,perms=out.permspikeavg,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr  pos ci.low ci.high   lod   pval
#K1RS.6110     1B.1  0.0      0     6.0 17.67 0.0000
#S2A_737386570   2A 73.5     68    77.0 16.57 0.0000
#S4B_547093481   4B 95.4     80   128.0  3.95 0.0174
#c4D.loc76       4D 76.0     69    80.8  8.82 0.0000
save(out.SIMLAavgspike, file = "D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-SIM-results")

# CIM
#window size 6; 4 covariate markers
out.CIMLAavgspike<-cim(SIMLAavgspike,pheno.col=8, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAavgspike,alpha=0.05,perms=out.permspikeavg,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod   pval
#c1A.loc75       1A  75.0   68.0    93.0  4.90 0.0020
#c1B.1.loc1    1B.1   1.0    0.0     4.0 26.84 0.0000
#S2A_737386570   2A  73.5   70.0    75.7 20.19 0.0000
#c3A.loc162      3A 162.0  136.0   165.0  4.45 0.0050
#S4B_603520056   4B 112.0   79.5   124.0  3.81 0.0234
#S4D_489467997   4D  80.8   77.0    80.8 10.75 0.0000
#c6B.loc185      6B 185.0  179.0   187.7  6.16 0.0006


#save QTL results
save(out.SIMLAavgspike, file = "D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-SIM-results")
save(out.CIMLAavgspike, file = "D:/wax/Field2021/QTL-results/R20-R21-LA-wax-spike-avg-CIM-results")
write.table(out.SIMLAavgspike,file="D:/wax/Field2021/QTL-results/LA_R20_waxleafB_SIM.csv", sep = ",")
write.table(out.CIMLAavgspike,file="D:/wax/Field2021/QTL-results/LA_R20_waxleafB_CIM.csv", sep = ",")



#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

SIMSPIKEM.1B.1.QTL<-makeqtl(LAmaprobR21,chr="1B.1",pos=0,qtl.name="K1RS.6110",what="prob")
SIMSPIKEM.1B.1.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=SIMSPIKEM.1B.1.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(SIMSPIKEM.1B.1.QTL.Model)

SIMSPIKEM.2A.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=73.5,qtl.name="S2A_737386570",what="prob")
SIMSPIKEM.2A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=SIMSPIKEM.2A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(SIMSPIKEM.2A.QTL.Model)

SIMSPIKEM.4B.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=95.4,qtl.name="S4B_547093481",what="prob")
SIMSPIKEM.4B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=SIMSPIKEM.4B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(SIMSPIKEM.4B.QTL.Model)

SIMSPIKEM.4D.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=76.0,qtl.name="c4D.loc76",what="prob")
SIMSPIKEM.4D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=SIMSPIKEM.4D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(SIMSPIKEM.4D.QTL.Model)

## Fit a multiple-QTL model 
MSPIKESIMQTLs<- makeqtl(LAmaprobR21,chr=c("1B.1","2A","4B","4D"),pos=c(0,73.5,95.4,76),
                    qtl.name=c("SIMSPIKEM.1B.1.QTL","SIMSPIKEM.2A.QTL", "SIMSPIKEM.4B.QTL", "SIMSPIKEM.4D.QTL"),
                    what="prob")
MSPIKESIMMQTL.model<- fitqtl(LAmaprobR21,pheno.col=8,qtl=MSPIKESIMQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(MSPIKESIMMQTL.model)


# b) out.CIM
## Fit single QTL model
CIMSPIKEM.1A.QTL<-makeqtl(LAmaprobR21,chr="1A",pos=75,qtl.name="c1A.loc75",what="prob")
CIMSPIKEM.1A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.1A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.1A.QTL.Model)

CIMSPIKEM.1B.1.QTL<-makeqtl(LAmaprobR21,chr="1B.1",pos=1,qtl.name="c1B.1.loc1",what="prob")
CIMSPIKEM.1B.1.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.1B.1.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.1B.1.QTL.Model)

CIMSPIKEM.2A.QTL<-makeqtl(LAmaprobR21,chr="2A",pos=73.5,qtl.name="S2A_737386570",what="prob")
CIMSPIKEM.2A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.2A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.2A.QTL.Model)

CIMSPIKEM.3A.QTL<-makeqtl(LAmaprobR21,chr="3A",pos=162,qtl.name="c3A.loc162",what="prob")
CIMSPIKEM.3A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.3A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.3A.QTL.Model)

CIMSPIKEM.4B.QTL<-makeqtl(LAmaprobR21,chr="4B",pos=112,qtl.name="S4B_603520056",what="prob")
CIMSPIKEM.4B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.4B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.4B.QTL.Model)

CIMSPIKEM.4D.QTL<-makeqtl(LAmaprobR21,chr="4D",pos=80.8,qtl.name="S4D_489467997",what="prob")
CIMSPIKEM.4D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.4D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.4D.QTL.Model)

CIMSPIKEM.6B.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=185,qtl.name="c6B.loc185",what="prob")
CIMSPIKEM.6B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=8,qtl=CIMSPIKEM.6B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMSPIKEM.6B.QTL.Model)


## Fit a multiple-QTL model 
SPIKEMCIMQTLs<- makeqtl(LAmaprobR21,chr=c("1A", "1B.1", "2A", "3A", "4B", "4D", "6B"),pos=c(75,1,73.5,162,112,80.8,185),
                       qtl.name=c("CIMSPIKEM.1A.QTL","CIMSPIKEM.1B.1.QTL", "CIMSPIKEM.2A.QTL", "CIMSPIKEM.3A.QTL", "CIMSPIKEM.4B.QTL", "CIMSPIKEM.4D.QTL", "CIMSPIKEM.6B.QTL"),
                       what="prob")
SPIKEMCIMQTL.model<- fitqtl(LAmaprobR21,pheno.col=8,qtl=SPIKEMCIMQTLs,formula=y~Q1*Q2*Q3*Q4*Q5*Q6*Q7,method="hk")
summary(SPIKEMCIMQTL.model)



#### RALEIGH 20 HEADING DATE (JULIAN) ####

#column 10

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR20HDJ<-scanone(LAmaprobR21,pheno.col=10,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR20HDJ, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.76
#5%  3.32

#save LOD result
save(out.permLAR20HDJ, file = "D:/wax/Field2021/QTL-results/R20-LA-HDJulian-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R20-LA-HDJulian-LODperm")
load("D:/wax/Field2021/QTL-results/R20-LA-HDJulian-SIM-results")
load("D:/wax/Field2021/QTL-results/R20-LA-HDJulian-CIM-results")

# SIM
SIMLAR20HD<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR20HD<-scanone(SIMLAR20HD,pheno.col=10,model="normal",method="em")
summary(out.SIMLAR20HD,alpha=0.05,perms=out.permLAR20HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#c2B.loc88       2B  88     85      93 16.04 0.000
#S5A_570227648   5A 237    230     258  3.83 0.016
#vrn.B1_AGS2K  5B.2   0      0      16  6.13 0.000
#c7D.loc15       7D  15     10      26  8.66 0.000

# CIM
#window size 6; 4 covariate markers
out.CIMLAR20HD<-cim(SIMLAR20HD,pheno.col=10, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR20HD,alpha=0.05,perms=out.permLAR20HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}


#save QTL results
save(out.SIMLAR20HD, file = "D:/wax/Field2021/QTL-results/R20-LA-HDJulian-SIM-results")
save(out.CIMLAR20HD, file = "D:/wax/Field2021/QTL-results/R20-LA-HDJulian-CIM-results")
write.table(out.SIMLAR20HD,file="D:/wax/Field2021/QTL-results/LA_R20_HDJulian_SIM.csv", sep = ",")
write.table(out.CIMLAR20HD,file="D:/wax/Field2021/QTL-results/LA_R20_HDJulian_CIM.csv", sep = ",")


#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

R20HD.2B.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=88,qtl.name="c2B.loc88",what="prob")
R20HD.2B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.2B.QTL.Model)

R20HD.5A.QTL<-makeqtl(LAmaprobR21,chr="5A",pos=237,qtl.name="S5A_570227648",what="prob")
R20HD.5A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.5A.QTL.Model)

R20HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
R20HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.5B.2.QTL.Model)

R20HD.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=15,qtl.name="c7D.loc15",what="prob")
R20HD.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.7D.QTL.Model)

## Fit a multiple-QTL model 
R20HDQTLs<- makeqtl(LAmaprobR21,chr=c("2B","5A","5B.2","7D"),pos=c(88,237,0,15),
                    qtl.name=c("R20HD.2B.QTL","R20HD.5A.QTL","R20HD.5B.2.QTL", "R20HD.7D.QTL"),
                    what="prob")
R20HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HDQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(R20HDMQTL.model)


# b) out.CIM
## Fit single QTL model
R20HD.2B.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=88,qtl.name="c2B.loc88",what="prob")
R20HD.2B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.2B.QTL.Model)

R20HD.5A.QTL<-makeqtl(LAmaprobR21,chr="5A",pos=237,qtl.name="S5A_570227648",what="prob")
R20HD.5A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=R20HD.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.5A.QTL.Model)

CIMR20HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=4,qtl.name="c5B.2.loc4",what="prob")
CIMR20HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=CIMR20HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR20HD.5B.2.QTL.Model)

CIMR20HD.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=16,qtl.name="c7D.loc16",what="prob")
CIMR20HD.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=10,qtl=CIMR20HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR20HD.7D.QTL.Model)


## Fit a multiple-QTL model 
CIMR20HDQTLs<- makeqtl(LAmaprobR21,chr=c("2B","5A","5B.2","7D"),pos=c(88,237,4,16),
                    qtl.name=c("R20HD.2B.QTL","R20HD.5A.QTL","CIMR20HD.5B.2.QTL", "CIMR20HD.7D.QTL"),
                    what="prob")
CIMR20HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=10,qtl=CIMR20HDQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(CIMR20HDMQTL.model)


#### KINSTON 20 HEADING DATE (JULIAN) ####

#column 11

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAK20HDJ<-scanone(LAmaprobR21,pheno.col=11,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAK20HDJ, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.68
#5%  3.28

#save LOD result
save(out.permLAK20HDJ, file = "D:/wax/Field2021/QTL-results/K20-LA-HDJulian-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/K20-LA-HDJulian-LODperm")
load("D:/wax/Field2021/QTL-results/K20-LA-HDJulian-SIM-results")
load("D:/wax/Field2021/QTL-results/K20-LA-HDJulian-CIM-results")

# SIM
SIMLAK20HD<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20HD<-scanone(SIMLAK20HD,pheno.col=11,model="normal",method="em")
summary(out.SIMLAK20HD,alpha=0.05,perms=out.permLAK20HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#vrn.B1_AGS2K 5B.2   0      0       4 33.08 0.000
#c6B.loc111     6B 111    105     127  4.98 0.003

# CIM
#window size 6; 4 covariate markers
out.CIMLAK20HD<-cim(SIMLAK20HD,pheno.col=11, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAK20HD,alpha=0.05,perms=out.permLAK20HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#S3B_53210923   3B  70   66.8      73  3.88 0.022
#c5A.loc239     5A 239  235.7     242 10.83 0.000
#vrn.B1_AGS2K 5B.2   0    0.0       4 47.29 0.000
#c7D.loc15      7D  15   14.0      21  8.84 0.000

#save QTL results
save(out.SIMLAK20HD, file = "D:/wax/Field2021/QTL-results/K20-LA-HDJulian-SIM-results")
save(out.CIMLAK20HD, file = "D:/wax/Field2021/QTL-results/K20-LA-HDJulian-CIM-results")
write.table(out.SIMLAK20HD,file="D:/wax/Field2021/QTL-results/LA_K20_HDJulian_SIM.csv", sep = ",")
write.table(out.CIMLAK20HD,file="D:/wax/Field2021/QTL-results/LA_K20_HDJulian_CIM.csv", sep = ",")



#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

K20HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
K20HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=K20HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K20HD.5B.2.QTL.Model)

K20HD.6B.QTL<-makeqtl(LAmaprobR21,chr="6B",pos=111,qtl.name="c6B.loc111",what="prob")
K20HD.6B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=K20HD.6B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K20HD.6B.QTL.Model)

## Fit a multiple-QTL model 
K20HDQTLs<- makeqtl(LAmaprobR21,chr=c("5B.2","6B"),pos=c(0,111),
                    qtl.name=c("K20HD.5B.2.QTL", "K20HD.6B.QTL"),
                    what="prob")
K20HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=11,qtl=K20HDQTLs,formula=y~Q1*Q2,method="hk")
summary(K20HDMQTL.model)


# b) out.CIM
## Fit single QTL model
CIMK20HD.3B.QTL<-makeqtl(LAmaprobR21,chr="3B",pos=70,qtl.name="S3B_53210923",what="prob")
CIMK20HD.3B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=CIMK20HD.3B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20HD.3B.QTL.Model)

CIMK20HD.5A.QTL<-makeqtl(LAmaprobR21,chr="5A",pos=239,qtl.name="c5A.loc239",what="prob")
CIMK20HD.5A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=CIMK20HD.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20HD.5A.QTL.Model)

CIMK20HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
CIMK20HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=CIMK20HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20HD.5B.2.QTL.Model)

CIMK20HD.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=15,qtl.name="c7D.loc15",what="prob")
CIMK20HD.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=11,qtl=CIMK20HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20HD.7D.QTL.Model)


## Fit a multiple-QTL model 
CIMK20HDQTLs<- makeqtl(LAmaprobR21,chr=c("3B","5A","5B.2","7D"),pos=c(70,239,0,15),
                       qtl.name=c("CIMK20HD.3B.QTL","CIMK20HD.5A.QTL","CIMK20HD.5B.2.QTL", "CIMK20HD.7D.QTL"),
                       what="prob")
CIMRK0HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=11,qtl=CIMK20HDQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(CIMRK0HDMQTL.model)



#### KINSTON 20 WINTER DORMANCY RELEASE (WDR) ####

#column 12

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAK20WDR<-scanone(LAmaprobR21,pheno.col=12,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAK20WDR, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.73
#5%  3.40

#save LOD result
save(out.permLAK20WDR, file = "D:/wax/Field2021/QTL-results/K20-LA-WDR-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/K20-LA-WDR-LODperm")
load("D:/wax/Field2021/QTL-results/K20-LA-WDR-SIM-results")
load("D:/wax/Field2021/QTL-results/K20-LA-WDR-CIM-results")

# SIM
SIMLAK20WDR<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20WDR<-scanone(SIMLAK20WDR,pheno.col=12,model="normal",method="em")
summary(out.SIMLAK20WDR,alpha=0.05,perms=out.permLAK20WDR,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high  lod pval
#vrn.B1_AGS2K 5B.2   0      0       4 34.9    0

# CIM
#window size 6; 4 covariate markers
out.CIMLAK20WDR<-cim(SIMLAK20WDR,pheno.col=12, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAK20WDR,alpha=0.05,perms=out.permLAK20WDR,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod  pval
#c3B.loc84       3B  84     81   88.00  3.87 0.013
#S5A_569871691   5A 236    234  240.00 11.18 0.000
#vrn.B1_AGS2K  5B.2   0      0    2.33 49.54 0.000
#c7D.loc15       7D  15     14   21.00  7.39 0.000

#save QTL results
save(out.SIMLAK20WDR, file = "D:/wax/Field2021/QTL-results/K20-LA-WDR-SIM-results")
save(out.CIMLAK20WDR, file = "D:/wax/Field2021/QTL-results/K20-LA-WDR-CIM-results")
write.table(out.SIMLAK20WDR,file="D:/wax/Field2021/QTL-results/LA_K20_WDR_SIM.csv", sep = ",")
write.table(out.CIMLAK20WDR,file="D:/wax/Field2021/QTL-results/LA_K20_WDR_CIM.csv", sep = ",")



#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

K20WDR.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
K20WDR.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=12,qtl=K20WDR.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K20WDR.5B.2.QTL.Model)


# b) out.CIM
## Fit single QTL model
CIMK20WDR.3B.QTL<-makeqtl(LAmaprobR21,chr="3B",pos=84,qtl.name="c3B.loc84",what="prob")
CIMK20WDR.3B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=12,qtl=CIMK20WDR.3B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20WDR.3B.QTL.Model)

CIMK20WDR.5A.QTL<-makeqtl(LAmaprobR21,chr="5A",pos=236,qtl.name="S5A_569871691",what="prob")
CIMK20WDR.5A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=12,qtl=CIMK20WDR.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20WDR.5A.QTL.Model)

CIMK20WDR.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
CIMK20WDR.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=12,qtl=CIMK20WDR.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20WDR.5B.2.QTL.Model)

CIMK20WDR.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=15,qtl.name="c7D.loc15",what="prob")
CIMK20WDR.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=12,qtl=CIMK20WDR.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMK20WDR.7D.QTL.Model)


## Fit a multiple-QTL model 
CIMK20WDRQTLs<- makeqtl(LAmaprobR21,chr=c("3B","5A","5B.2","7D"),pos=c(84,236,0,15),
                       qtl.name=c("CIMK20WDR.3B.QTL","CIMK20WDR.5A.QTL","CIMK20WDR.5B.2.QTL", "CIMK20WDR.7D.QTL"),
                       what="prob")
CIMRK20WDRMQTL.model<- fitqtl(LAmaprobR21,pheno.col=12,qtl=CIMK20WDRQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(CIMRK20WDRMQTL.model)

#### RALEIGH 21 HEADING DATE (JULIAN) REDO ####

#column 4

#determine LOD threshold
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR21HDJ<-scanone(LAmaprobR21,pheno.col=4,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR21HDJ, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.73
#5%  3.36

#save LOD result
save(out.permLAR21HDJ, file = "D:/wax/Field2021/QTL-results/R21-LA-HDJulian-REDO-LODperm")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/wax/Field2021/QTL-results/R21-LA-HDJulian-REDO-LODperm")
load("D:/wax/Field2021/QTL-results/R21-LA-HDJulian-SIM-results")
load("D:/wax/Field2021/QTL-results/R21-LA-HDJulian-CIM-results")

# SIM
SIMLAR21HDJ<-calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR21HDJ<-scanone(SIMLAR21HDJ,pheno.col=4,model="normal",method="em")
summary(out.SIMLAR21HDJ,alpha=0.05,perms=out.permLAR21HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high   lod pval
#c2B.loc92      2B  92     85      95  8.05    0
#vrn.B1_AGS2K 5B.2   0      0       4 16.78    0
#c7D.loc13      7D  13      4      26  7.93    0

# CIM
#window size 6; 4 covariate markers
out.CIMLAR21HDJ<-cim(SIMLAR21HDJ,pheno.col=4, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAR21HDJ,alpha=0.05,perms=out.permLAR21HDJ,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod  pval
#c2B.loc91       2B  91.0     89    94.0  8.04 0.000
#S5A_570528220   5A 238.0    237   244.0  6.29 0.001
#vrn.B1_AGS2K  5B.2   0.0      0     4.0 28.52 0.000
#S6A_541404429   6A  75.3     68    81.6  4.27 0.010
#c7D.loc15       7D  15.0     14    21.0 17.54 0.000

#save QTL results
save(out.SIMLAR21HDJ, file = "D:/wax/Field2021/QTL-results/R21-LA-HDJulian-SIM-results")
save(out.CIMLAR21HDJ, file = "D:/wax/Field2021/QTL-results/R21-LA-HDJulian-CIM-results")
write.table(out.SIMLAR21HDJ,file="D:/wax/Field2021/QTL-results/LA_R21_HDJulian_SIM.csv", sep = ",")
write.table(out.CIMLAR21HDJ,file="D:/wax/Field2021/QTL-results/LA_R21_HDJulian_CIM.csv", sep = ",")



#### MODEL PVE AND EPISTASIS 

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

R21HD.2B.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=92,qtl.name="c2B.loc92",what="prob")
R21HD.2B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HD.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21HD.2B.QTL.Model)

R21HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
R21HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21HD.5B.2.QTL.Model)

R21HD.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=13,qtl.name="c7D.loc13",what="prob")
R21HD.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21HD.7D.QTL.Model)

## Fit a multiple-QTL model 
R21HDQTLs<- makeqtl(LAmaprobR21,chr=c("2B","5B.2","7D"),pos=c(92,0,13),
                       qtl.name=c("R21HD.2B.QTL","R21HD.5B.2.QTL", "R21HD.7D.QTL"),
                       what="prob")
R21HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HDQTLs,formula=y~Q1*Q2*Q3,method="hk")
summary(R21HDMQTL.model)

# b) out.CIM
## Fit single QTL model
CIMR21HD.2B.QTL<-makeqtl(LAmaprobR21,chr="2B",pos=91,qtl.name="c2B.loc91",what="prob")
CIMR21HD.2B.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=CIMR21HD.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR21HD.2B.QTL.Model)

CIMR21HD.5A.QTL<-makeqtl(LAmaprobR21,chr="5A",pos=238,qtl.name="S5A_570528220",what="prob")
CIMR21HD.5A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=CIMR21HD.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR21HD.5A.QTL.Model)

R21HD.5B.2.QTL<-makeqtl(LAmaprobR21,chr="5B.2",pos=0,qtl.name="vrn.B1_AGS2K",what="prob")
R21HD.5B.2.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21HD.5B.2.QTL.Model)

CIMR21HD.6A.QTL<-makeqtl(LAmaprobR21,chr="6A",pos=75.3,qtl.name="S6A_541404429",what="prob")
CIMR21HD.6A.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=CIMR21HD.6A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR21HD.6A.QTL.Model)

CIMR21HD.7D.QTL<-makeqtl(LAmaprobR21,chr="7D",pos=15,qtl.name="c7D.loc15",what="prob")
CIMR21HD.7D.QTL.Model<-fitqtl(LAmaprobR21,pheno.col=4,qtl=CIMR21HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR21HD.7D.QTL.Model)


## Fit a multiple-QTL model 
R21HDCIMQTLs<- makeqtl(LAmaprobR21,chr=c("2B", "5A", "5B.2", "6A", "7D"),pos=c(91,238,0,75.3,15),
                          qtl.name=c("CIMR21HD.2B.QTL","CIMR21HD.5A.QTL", "R21HD.5B.2.QTL", "CIMR21HD.6A.QTL", "CIMR21HD.7D.QTL"),
                          what="prob")
R21HDMCIMQTL.model<- fitqtl(LAmaprobR21,pheno.col=4,qtl=R21HDCIMQTLs,formula=y~Q1*Q2*Q3*Q4*Q5,method="hk")
summary(R21HDMCIMQTL.model)
