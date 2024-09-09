# Hilliard x GA DH pop - QTL Mapping - Raleigh+Kinston 2021 - Wax/Glaucousness Phenotypes - Leaf + Heads - Heading Date

#options(max.print=999999)
#memory.limit(size = NA)
library("qtl")
library("magrittr")
library("ASMap")
#library("dplyr")

##### MERGE PHENOTYPE AND GENOTYPE DATA #####

library(dplyr)
setwd("C:/Users/Daniela/Documents/Wax/HilliardxGA/")
pheno <- read.csv("HG_R21_K21_waxtrts.csv")
head(pheno)

geno <- read.csv("UX1992-HilliardxGA-RIL-Miller-linkage-map.csv")
geno[1:10,1:10]

joined_df <- merge(pheno, geno, by.x = "genotype", 
                   by.y = "genotype", all.x = TRUE, all.y = TRUE)
dim(pheno)
dim(geno)
dim(joined_df)

setwd("~/Wax/HGfinal/")
write.csv(joined_df, "HilliardxGA-RIL-linkage-map_waxtrts_Miller.csv")


####  READ IN MAP FILE ####

map<-read.cross(format="csv",file="~/Wax/HGfinal/HilliardxGA-RIL-linkage-map_waxtrts_Miller.csv",estimate.map=FALSE,
                  na.strings=c(".","-","NA"),genotypes=c("AA","BB"),crosstype="riself")%>%jittermap(.)%>%calc.genoprob(map.function="kosambi")
summary(map)


#### ADD 3A MARKERS TO MAP ####

## read in 3A markers
KASP<-read.cross(format="csv",file="HG_3A_markers_format.csv",estimate.map=FALSE,
                 na.strings=c("NA"),genotypes=c("A","B"),
                 crosstype="riself")%T>% {print(summary(.))}
summary(KASP)


together<-combineMap(map,KASP,id="genotype",keep.all=TRUE,
                        merge.by = "genotype")
jittermap(together)
summary(together)
class(together)

#-Build the map
geneticMap<-mstmap(together,pop.type="ARIL",dist.fun="kosambi",
                           id="genotype",objective.fun="COUNT",noMap.dist=15,noMap.size=0,
                           bychr=TRUE,miss.thresh=1,mvest.bc=FALSE,detectBadData=FALSE, 
                           as.cross=TRUE, return.imputed=FALSE,trace=TRUE,anchor=TRUE,
                           p.value=1)


nmar(geneticMap)
chrlen(geneticMap)
summaryMap<-summary.map(geneticMap)%T>%{print(.)}

#final map used for QTL mapping:
write.cross(geneticMap,"csv",filestem="UX1992-HilliardxGA-RIL-Miller-3Alinkage-map",
            c("1A","1B","1D","2A","2B","2D","3A","3B",
              "4A","4B","5A","5B","6A","6B","6D",
              "7A","7B.1","7B.2","7D"))

#### SPIKE WAX - COMBINED ENV MODEL ####

#load final map from file
setwd("C:/Users/Daniela/Documents/Wax/HGfinal/")
#map <- geneticMap
map<-read.cross(format="csv",file="UX1992-HilliardxGA-RIL-Miller-3Alinkage-map.csv",estimate.map=FALSE,
                na.strings=c(".","-","NA"),genotypes=c("A","B"),crosstype="riself")%>%jittermap(.)%>%calc.genoprob(map.function="kosambi")
summary(map)


#column 9 (spike_emmeans)

##determine LOD threshold
maprob<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permspike<-scanone(maprob,pheno.col=9,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permspike, alpha=c(0.20,0.05))

##LOD thresholds (1000 permutations)
##lod
#20% 2.49
#5%  3.17

##save LOD result
setwd("~/Wax/HGfinal/QTL-results")
save(out.permspike, file = "HG-spikewaxcomb-LODperm1000")

#after first time, quick load:
maprob<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("HG-spikewaxcomb-LODperm1000")
load("HG-spikewaxcomb-SIM-results")
load("HG-spikewaxcomb-CIM-results")

# SIM
SIMspike<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMspike<-scanone(SIMspike,pheno.col=9,model="normal",method="em")
summary(out.SIMspike,alpha=0.05,perms=out.permspike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#             chr   pos ci.low ci.high   lod  pval
#S1A_6370170   1A   0.0      0      15  3.69 0.013
#c3A.loc119    3A 119.0    111     129 10.65 0.000
#S5B_81171755  5B  39.7     27      51  3.90 0.010
#c7A.loc192    7A 192.0    175     207  3.35 0.035


# CIM
#window size 6; 4 covariate markers
out.CIMspike<-cim(SIMspike,pheno.col=9, n.marcovar=4,method="em",window=6)
summary(out.CIMspike,alpha=0.05,perms=out.permspike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#               chr   pos ci.low ci.high   lod  pval
#S1A_6370170    1A   0.0      0       4  4.95 0.001
#S3A_621010472  3A 126.5    122     129 12.89 0.000
#S5B_78423903   5B  37.9     34      41  5.04 0.001
#S7A_13777429   7A  16.8      4     193  3.19 0.049


#save QTL results
save(out.SIMspike, file = "HG-spikewaxcomb-SIM-results")
save(out.CIMspike, file = "HG-spikewaxcomb-CIM-results")
write.table(out.SIMspike,file="HG-spikewaxcomb_SIM.csv", sep = ",")
write.table(out.CIMspike,file="HG-spikewaxcomb_CIM.csv", sep = ",")

#### SPIKE COMBINED - MODEL PVE AND EPISTASIS ####


# a) SIM - Fit single QTL model

#1A
c1A.SIM.spike.QTL<-makeqtl(maprob,chr="1A",pos=0,qtl.name="S1A_6370170",what="prob")
c1A.SIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c1A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c1A.SIM.spike.QTL.Model)

#3A
c3A.SIM.spike.QTL<-makeqtl(maprob,chr="3A",pos=119,qtl.name="c3A.loc119",what="prob")
c3A.SIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c3A.SIM.spike.QTL.Model)

#5B
c5B.SIM.spike.QTL<-makeqtl(maprob,chr="5B",pos=39.7,qtl.name="S5B_81171755",what="prob")
c5B.SIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c5B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c5B.SIM.spike.QTL.Model)

#7A
c7A.SIM.spike.QTL<-makeqtl(maprob,chr="7A",pos=192,qtl.name="c7A.loc192",what="prob")
c7A.SIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c7A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c7A.SIM.spike.QTL.Model)

#three QTL changed position in CIM:

#3A
c3A.CIM.spike.QTL<-makeqtl(maprob,chr="3A",pos=126.5,qtl.name="S3A_621010472",what="prob")
c3A.CIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c3A.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c3A.CIM.spike.QTL.Model)

#5B
c5B.CIM.spike.QTL<-makeqtl(maprob,chr="5B",pos=37.9,qtl.name="S5B_78423903",what="prob")
c5B.CIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c5B.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c5B.CIM.spike.QTL.Model)

#7A
c7A.CIM.spike.QTL<-makeqtl(maprob,chr="7A",pos=16.8,qtl.name="S7A_13777429",what="prob")
c7A.CIM.spike.QTL.Model<-fitqtl(maprob,pheno.col=9,qtl=c7A.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c7A.CIM.spike.QTL.Model)


#CIM
## Fit a multiple-QTL model 
spikeQTL<- makeqtl(maprob,chr=c("1A","3A","5B","7A"),pos=c(0.0,126.5,37.9,16.8),
                    qtl.name=c("1A.spike.QTL","3A.spike.QTL","5B.spike.QTL","7A.spike.QTL"),
                    what="prob")
spikeQTLs.model<- fitqtl(maprob,pheno.col=9,qtl=spikeQTL,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(spikeQTLs.model)




# OLD - SIM

## Fit a multiple-QTL model 
spikeQTLs<- makeqtl(maprob,chr=c("1A","3A","5B","7A"),pos=c(0,119,39.7,192),
                       qtl.name=c("c1A.SIM.spike.QTL","c3A.SIM.spike.QTL","c5B.SIM.spike.QTL","c7A.SIM.spike.QTL"),
                       what="prob")
spikeQTL.model<- fitqtl(maprob,pheno.col=9,qtl=spikeQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(spikeQTL.model)

#### No Significant epistasis 


#### LEAF WAX COMBINED ####

#column 8 (leaf_emmeans)

# leafmeans ##determine LOD threshold
maprob<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permleaf<-scanone(maprob,pheno.col=8,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permleaf, alpha=c(0.20,0.05))
#LOD thresholds (100 permutations)
#20% 2.50
#5%  3.13

##save LOD result
setwd("~/Wax/HGfinal/QTL-results")
save(out.permspike, file = "HG-leafwaxcomb-LODperm1000")

#after first time, quick load:
maprob<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("HG-leafwaxcomb-LODperm1000")
load("HG-leafwaxcomb-SIM-results")
load("HG-leafwaxcomb-CIM-results")

# SIM
SIMleaf<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMleaf<-scanone(SIMleaf,pheno.col=8,model="normal",method="em")
summary(out.SIMleaf,alpha=0.05,perms=out.permleaf,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr   pos ci.low ci.high   lod  pval
#S2D_17164152   2D 2.00e-06      0       6  4.64 0.002
#c3A.loc119     3A 1.19e+02    112     125 17.73 0.000
#S4B_559855530  4B 1.50e+01      2      30  3.44 0.024

# CIM
#window size 6; 4 covariate markers
out.CIMleaf<-cim(SIMleaf,pheno.col=8, n.marcovar=4,method="em",window=6)
summary(out.CIMleaf,alpha=0.05,perms=out.permleaf,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#               chr   pos ci.low ci.high   lod  pval
#c1A.loc84      1A 8.40e+01   78.4    92.4  3.66 0.011
#S2A_1872142    2A 0.00e+00    0.0     4.0  3.63 0.011
#S2D_17164152   2D 2.00e-06    0.0     3.8  6.98 0.000
#c3A.loc122     3A 1.22e+02  120.9   126.0 22.54 0.000
#S4B_559855530  4B 1.50e+01   12.0    18.8  3.97 0.006


#save QTL results
save(out.SIMleaf, file = "HG-leafwaxcomb-SIM-results")
save(out.CIMleaf, file = "HG-leafwaxcomb-CIM-results")
write.table(out.SIMleaf,file="HG-leafewaxcomb_SIM.csv", sep = ",")
write.table(out.CIMleaf,file="HG-leafwaxcomb_CIM.csv", sep = ",")


#### LEAF COMBINED - MODEL PVE AND EPISTASIS ####

# a) SIM - Fit single QTL model

#2D
c2D.SIM.leaf.QTL<-makeqtl(maprob,chr="2D",pos=2.00e-06,qtl.name="S2D_17164152",what="prob")
c2D.SIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c2D.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c2D.SIM.leaf.QTL.Model)

#3A
c3A.SIM.leaf.QTL<-makeqtl(maprob,chr="3A",pos=119,qtl.name="c3A.loc119",what="prob")
c3A.SIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c3A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c3A.SIM.leaf.QTL.Model)

#4B
c4B.SIM.leaf.QTL<-makeqtl(maprob,chr="4B",pos=1.50e+01,qtl.name="S4B_559855530",what="prob")
c4B.SIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c4B.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c4B.SIM.leaf.QTL.Model)


#three QTL changed position in CIM:

#1A
c1A.CIM.leaf.QTL<-makeqtl(maprob,chr="1A",pos=84,qtl.name="c1A.loc84",what="prob")
c1A.CIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c1A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c1A.CIM.leaf.QTL.Model)

#2A
c2A.CIM.leaf.QTL<-makeqtl(maprob,chr="2A",pos=0,qtl.name="S2A_1872142",what="prob")
c2A.CIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c2A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c2A.CIM.leaf.QTL.Model)

#2D
c2D.CIM.leaf.QTL<-makeqtl(maprob,chr="2D",pos=2.00e-06,qtl.name="S2D_17164152",what="prob")
c2D.CIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c2D.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c2D.CIM.leaf.QTL.Model)

#3A
c3A.CIM.leaf.QTL<-makeqtl(maprob,chr="3A",pos=122,qtl.name="c3A.loc122",what="prob")
c3A.CIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c3A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c3A.CIM.leaf.QTL.Model)

#4B
c4B.CIM.leaf.QTL<-makeqtl(maprob,chr="4B",pos=1.50e+01,qtl.name="S4B_559855530",what="prob")
c4B.CIM.leaf.QTL.Model<-fitqtl(maprob,pheno.col=8,qtl=c4B.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(c4B.CIM.leaf.QTL.Model)


#CIM
## Fit a multiple-QTL model 
leafQTL<- makeqtl(maprob,chr=c("1A","2A","2D","3A","4B"),pos=c(84,0,2.00e-06,122,1.50e+01),
                   qtl.name=c("1A.leaf.QTL","2A.leaf.QTL","2D.leaf.QTL","3A.leaf.QTL","4B.leaf.QTL"),
                   what="prob")
leafQTLs.model<- fitqtl(maprob,pheno.col=8,qtl=leafQTL,formula=y~Q1*Q2*Q3*Q4*Q5,method="hk")
summary(leafQTLs.model)


# OLD - SIM

## Fit a multiple-QTL model 
leafQTLs<- makeqtl(maprob,chr=c("2D","3A","4B"),pos=c(2.00e-06,119,1.50e+01),
                    qtl.name=c("c2D.SIM.leaf.QTL","c3A.SIM.leaf.QTL","c4B.SIM.leaf.QTL"),
                    what="prob")
leafQTL.model<- fitqtl(maprob,pheno.col=8,qtl=leafQTLs,formula=y~Q1*Q2*Q3,method="hk")
summary(leafQTL.model)

#### No Significant epistasis 

#### RALEIGH 21 SPIKE WAX ####

#map <- geneticMap
map<-read.cross(format="csv",file="UX1992-HilliardxGA-RIL-Miller-3Alinkage-map.csv",estimate.map=FALSE,
                na.strings=c(".","-","NA"),genotypes=c("A","B"),crosstype="riself")%>%jittermap(.)%>%calc.genoprob(map.function="kosambi")
summary(map)


#column 2 (first true phenotype)

##determine LOD threshold
maprobR21<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
out.permLAR21spike<-scanone(maprobR21,pheno.col=2,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permLAR21spike, alpha=c(0.20,0.05))

##LOD thresholds (1000 permutations)
##lod
#20% 2.58
#5%  3.22

##save LOD result
save(out.permLAR21spike, file = "R21-HG-spikewax-LODperm1000")

#after first time, quick load:
maprobR21<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("R21-HG-spikewax-LODperm1000")
load("R21-spikewax-SIM-results")
load("R21-spikewax-CIM-results")

# SIM
SIMLAK20WDR<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20WDR<-scanone(SIMLAK20WDR,pheno.col=2,model="normal",method="em")
summary(out.SIMLAK20WDR,alpha=0.05,perms=out.permLAR21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#             chr   pos ci.low ci.high   lod  pval
#S1A_6370170  1A   0      0      12  3.25 0.044
#c3A.loc114   3A 114    111     127 15.16 0.000
#c5B.loc43    5B  43     25      57  3.58 0.021


# CIM
#window size 6; 4 covariate markers
out.CIMLAK20WDR<-cim(SIMLAK20WDR,pheno.col=2, n.marcovar=4,method="em",window=6)
summary(out.CIMLAK20WDR,alpha=0.05,perms=out.permLAR21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#               chr   pos ci.low ci.high   lod  pval
#S3A_617060163  3A 125.0    120     127 17.01    0
#S5B_81171755   5B  39.7     36      43  5.13    0


#save QTL results
save(out.SIMLAK20WDR, file = "R21-spikewax-SIM-results")
save(out.CIMLAK20WDR, file = "R21-spikewax-CIM-results")
write.table(out.SIMLAK20WDR,file="R21-spikewax_SIM.csv", sep = ",")
write.table(out.CIMLAK20WDR,file="R21-spikewax_CIM.csv", sep = ",")

#plot
#significant QTL LGs only 
plot(out.SIMLAK20WDR,out.CIMLAK20WDR, chr=c("1A", "3A", "5B"), col=c("cyan", "blue"), main="HilliardxGA06493-13LE6 RIL Spike Glaucousness Raleigh 2021")
add.threshold(out.SIMLAK20WDR,perms=out.permLAR21spike,alpha=0.05, col="orange")
plot(out.SIMLAK20WDR,out.CIMLAK20WDR, chr=c("1A", "3A", "5B"), col=c("cyan", "blue"), main="Spike Glaucousness Raleigh 2021")
add.threshold(out.SIMLAK20WDR,perms=out.permLAR21spike,alpha=0.05, col="orange")


#### R21 SPIKE WAX - MODEL PVE AND EPISTASIS ####


#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#1A
R21.1A.SIM.spike.QTL<-makeqtl(maprobR21,chr="1A",pos=0,qtl.name="S1A_6370170",what="prob")
R21.1A.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=2,qtl=R21.1A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1A.SIM.spike.QTL.Model)

#3A
R21.3A.SIM.spike.QTL<-makeqtl(maprobR21,chr="3A",pos=114,qtl.name="c3A.loc114",what="prob")
R21.3A.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=2,qtl=R21.3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.SIM.spike.QTL.Model)

#5B
R21.5B.SIM.spike.QTL<-makeqtl(maprobR21,chr="5B",pos=43,qtl.name="c5B.loc43",what="prob")
R21.5B.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=2,qtl=R21.5B.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.5B.SIM.spike.QTL.Model)


#only one QTL changed position in CIM:

#3A
R21.3A.CIM.spike.QTL<-makeqtl(maprobR21,chr="3A",pos=125,qtl.name="S3A_617060163",what="prob")
R21.3A.CIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=2,qtl=R21.3A.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.CIM.spike.QTL.Model)

#5B
R21.5B.CIM.spike.QTL<-makeqtl(maprobR21,chr="5B",pos=39.7,qtl.name="S5B_81171755",what="prob")
R21.5B.CIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=2,qtl=R21.5B.CIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.5B.CIM.spike.QTL.Model)

## Fit a multiple-QTL model 
R21spikeQTLs<- makeqtl(maprobR21,chr=c("1A","3A","5B"),pos=c(0,114,43),
                        qtl.name=c("R21.1A.SIM.spike.QTL","R21.3A.SIM.spike.QTL","R21.5B.SIM.spike.QTL"),
                        what="prob")
R21spikeQTL.model<- fitqtl(maprobR21,pheno.col=2,qtl=R21spikeQTLs,formula=y~Q1*Q2*Q3,method="hk")
summary(R21spikeQTL.model)

#### No Significant epistasis 



#### KINSTON 21 SPIKE WAX ####

#column 3 (SECOND true phenotype)

#determine LOD threshold
out.permK21spike<-scanone(maprobR21,pheno.col=3,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permK21spike, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.53
#5%  3.12

#save LOD result
save(out.permK21spike, file = "K21-HG-spikewax-LODperm1000")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("K21-HG-spikewax-LODperm1000")
load("K21-spikewax-SIM-results")
load("K21-spikewax-CIM-results")


# SIM
SIMLAK20WDR<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20WDR<-scanone(SIMLAK20WDR,pheno.col=3,model="normal",method="em")
summary(out.SIMLAK20WDR,alpha=0.05,perms=out.permK21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high  lod  pval
#c3A.loc119  3A 119     64     143 4.45 0.002



# CIM
#window size 6; 4 covariate markers
out.CIMLAK20WDR<-cim(SIMLAK20WDR,pheno.col=3, n.marcovar=4,method="em",window=6)
summary(out.CIMLAK20WDR,alpha=0.05,perms=out.permK21spike,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high  lod  pval
#c3A.loc118  3A 118    113     120 5.07 0.001



#save QTL results
save(out.SIMLAK20WDR, file = "K21-spikewax-SIM-results")
save(out.CIMLAK20WDR, file = "K21-spikewax-CIM-results")
write.table(out.SIMLAK20WDR,file="K21-spikewax_SIM.csv", sep = ",")
write.table(out.CIMLAK20WDR,file="K21-spikewax_CIM.csv", sep = ",")


#plot
#significant QTL LGs only 
plot(out.SIMLAK20WDR,out.CIMLAK20WDR, chr=c("1A", "3A"), col=c("cyan", "blue"), main="Spike Glaucousness Kinston 2021")
add.threshold(out.SIMLAK20WDR,perms=out.permK21spike,alpha=0.05, col="orange")


#### K21 SPIKE WAX - MODEL PVE AND EPISTASIS  ####

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#1A
K21.1A.SIM.spike.QTL<-makeqtl(maprobR21,chr="1A",pos=18.8,qtl.name="S1A_9236441",what="prob")
K21.1A.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=3,qtl=K21.1A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.1A.SIM.spike.QTL.Model)

#3A
K21.3A.SIM.spike.QTL<-makeqtl(maprobR21,chr="3A",pos=119,qtl.name="c3A.loc119",what="prob")
K21.3A.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=3,qtl=K21.3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.3A.SIM.spike.QTL.Model)

#3A CIM
K21.3A.SIM.spike.QTL<-makeqtl(maprobR21,chr="3A",pos=118,qtl.name="c3A.loc118",what="prob")
K21.3A.SIM.spike.QTL.Model<-fitqtl(maprobR21,pheno.col=3,qtl=K21.3A.SIM.spike.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.3A.SIM.spike.QTL.Model)



## Fit a multiple-QTL model 
K21spikeQTLs<- makeqtl(maprobR21,chr=c("1A","3A"),pos=c(18.8,69),
                       qtl.name=c("K21.1A.SIM.spike.QTL","K21.3A.SIM.spike.QTL"),
                       what="prob")
K21spikeQTL.model<- fitqtl(maprobR21,pheno.col=3,qtl=K21spikeQTLs,formula=y~Q1*Q2,method="hk")
summary(K21spikeQTL.model)

#### No Significant epistasis 



#### RALEIGH 21 LEAF BOTTOM WAX ####

#column 4

#determine LOD threshold
out.permR21leafB<-scanone(maprobR21,pheno.col=4,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permR21leafB, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.49
#5%  3.16

#save LOD result
save(out.permR21leafB, file = "R21-HG-leafBwax-LODperm1000")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("R21-HG-leafBwax-LODperm1000")
load("R21-leafwax-SIM-results")
load("R21-leafwax-CIM-results")



# SIM
SIMLAR20leafB<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAR20leafB<-scanone(SIMLAR20leafB,pheno.col=4,model="normal",method="em")
summary(out.SIMLAR20leafB,alpha=0.05,perms=out.permR21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr pos ci.low ci.high  lod pval
#c3A.loc116  3A 116    112     124 18.3    0


# CIM
#window size 6; 4 covariate markers
out.CIMLAR20leafB<-cim(SIMLAR20leafB,pheno.col=4, n.marcovar=4,method="em",window=6)
summary(out.CIMLAR20leafB,alpha=0.05,perms=out.permR21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr    pos ci.low ci.high   lod  pval
#c2A.loc3       2A   3      1     8.0  4.22 0.007
#c3A.loc116     3A 116    113   120.0 22.37 0.000
#S4B_559855530  4B  15     12    18.8  3.66 0.022


#save QTL results
save(out.SIMLAR20leafB, file = "R21-leafwax-SIM-results")
save(out.CIMLAR20leafB, file = "R21-leafwax-CIM-results")
write.table(out.SIMLAR20leafB,file="R21-leafwax_SIM.csv", sep = ",")
write.table(out.CIMLAR20leafB,file="R21-leafwax_SIM.csv", sep = ",")


#### 4: Plots

#SIM plots
plot(out.SIMLAR20leafB, main = "Standard Interval Mapping - LA R21 wax leaf bottom")
add.threshold(out.SIMLAR20leafB,perms=out.permR21leafB,alpha=0.05, col="green")

#CIM plots
plot(out.CIMLAR20leafB, main = "Composite Interval Mapping - LA R21 wax leaf bottom")
add.threshold(out.CIMLAR20leafB,perms=out.permR21leafB,alpha=0.05, col="blue")

#significant QTL LGs only 
plot(out.SIMLAR20leafB,out.CIMLAR20leafB, chr=c("1A", "3A", "4B"), col=c("cyan", "blue"), main="HilliardxGA06493-13LE6 RIL Leaf Glaucousness Raleigh 2021")
add.threshold(out.SIMLAR20leafB,perms=out.permR21leafB,alpha=0.05, col="orange")


#### R21 LEAF WAX - MODEL PVE AND EPISTASIS  ####

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

#3A
R21.3A.SIM.leaf.QTL<-makeqtl(maprobR21,chr="3A",pos=116,qtl.name="c3A.loc116",what="prob")
R21.3A.SIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=4,qtl=R21.3A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.SIM.leaf.QTL.Model)

# CIM models

#1A
R21.1A.SIM.leaf.QTL<-makeqtl(maprobR21,chr="1A",pos=87,qtl.name="S1A_520868230",what="prob")
R21.1A.SIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=4,qtl=R21.1A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.1A.SIM.leaf.QTL.Model)

#2A
R21.2A.CIM.leaf.QTL<-makeqtl(maprobR21,chr="2A",pos=3,qtl.name="c2A.loc3",what="prob")
R21.2A.CIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=4,qtl=R21.2A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.2A.CIM.leaf.QTL.Model)

#3A
R21.3A.CIM.leaf.QTL<-makeqtl(maprobR21,chr="3A",pos=116,qtl.name="c3A.loc116",what="prob")
R21.3A.CIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=4,qtl=R21.3A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.3A.CIM.leaf.QTL.Model)

#4B
R21.4B.CIM.leaf.QTL<-makeqtl(maprobR21,chr="4B",pos=15,qtl.name="S4B_559855530",what="prob")
R21.4B.CIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=4,qtl=R21.4B.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21.4B.CIM.leaf.QTL.Model)




## Fit a multiple-QTL model 
R21leafQTLs<- makeqtl(maprobR21,chr=c("1A","2A","3A","4B"),pos=c(87,4.65,114,19.02),
                       qtl.name=c("R21.1A.SIM.leaf.QTL","R21.2A.CIM.leaf.QTL","R21.3A.CIM.leaf.QTL", "R21.4B.CIM.leaf.QTL"),
                       what="prob")
R21leafQTL.model<- fitqtl(maprobR21,pheno.col=4,qtl=R21leafQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(R21leafQTL.model)

#### No Significant epistasis 



#### KINSTON 21 LEAF BOTTOM WAX ####

#column 5 

#determine LOD threshold
out.permK21leafB<-scanone(maprobR21,pheno.col=5,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permK21leafB, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.56
#5%  3.22

#save LOD result
save(out.permK21leafB, file = "K21-HG-leafBwax-LODperm1000")

#after first time, quick load:
load("K21-HG-leafBwax-LODperm1000")
load("K21-leafwax-SIM-results")
load("K21-leafwax-CIM-results")

# SIM
SIMLAK20WDR<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20WDR<-scanone(SIMLAK20WDR,pheno.col=5,model="normal",method="em")
summary(out.SIMLAK20WDR,alpha=0.05,perms=out.permK21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr      pos ci.low ci.high   lod  pval
#c1B.loc85     1B  85     72    98.0  3.89 0.009
#S2D_17457432  2D   0      0     3.8  5.59 0.001
#c3A.loc122    3A 122    117   125.0 11.75 0.000

# CIM
#window size 6; 4 covariate markers
out.CIMLAK20WDR<-cim(SIMLAK20WDR,pheno.col=5, n.marcovar=4,method="em",window=6)#,map.function="kosambi",n.perm=1000)
summary(out.CIMLAK20WDR,alpha=0.05,perms=out.permK21leafB,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr      pos ci.low ci.high   lod  pval
#c1B.loc84     1B  84     83    90.0  4.84 0.001
#S2D_17457432  2D   0      0     3.8  4.73 0.002
#c3A.loc122    3A 122    118   123.3 17.05 0.000


#save QTL results
save(out.SIMLAK20WDR, file = "K21-leafwax-SIM-results")
save(out.CIMLAK20WDR, file = "K21-leafwax-CIM-results")
write.table(out.SIMLAK20WDR,file="K21-leafwax_SIM.csv", sep = ",")
write.table(out.CIMLAK20WDR,file="K21-leafwax_CIM.csv", sep = ",")

#plot
#significant QTL LGs only 
plot(out.SIMLAK20WDR,out.CIMLAK20WDR, chr=c("1B", "2D", "3A"), col=c("cyan", "blue"), main="Leaf Glaucousness Kinston 2021")
add.threshold(out.SIMLAK20WDR,perms=out.permR21leafB,alpha=0.05, col="orange")


#### K21 LEAF BOTTOM MODEL PVE AND EPISTASIS ####

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model

# a) SIM - Fit single QTL model

#1B
K21.1B.SIM.leaf.QTL<-makeqtl(maprobR21,chr="1B",pos=85,qtl.name="c1B.loc85",what="prob")
K21.1B.SIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=5,qtl=K21.1B.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.1B.SIM.leaf.QTL.Model)

#2D
K21.2D.SIM.leaf.QTL<-makeqtl(maprobR21,chr="2D",pos=2,qtl.name="S2D_17457432",what="prob")
K21.2D.SIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=5,qtl=K21.2D.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.2D.SIM.leaf.QTL.Model)

#3A
K21.3A.SIM.leaf.QTL<-makeqtl(maprobR21,chr="3A",pos=122,qtl.name="c3A.loc122",what="prob")
K21.3A.SIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=5,qtl=K21.3A.SIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.3A.SIM.leaf.QTL.Model)

# CIM models

#1B
K21.1B.CIM.leaf.QTL<-makeqtl(maprobR21,chr="1B",pos=84,qtl.name="c1B.loc84",what="prob")
K21.1B.CIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=5,qtl=K21.1B.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.1B.CIM.leaf.QTL.Model)

#3A
K21.3A.CIM.leaf.QTL<-makeqtl(maprobR21,chr="3A",pos=114,qtl.name="c3A.loc114",what="prob")
K21.3A.CIM.leaf.QTL.Model<-fitqtl(maprobR21,pheno.col=5,qtl=K21.3A.CIM.leaf.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(K21.3A.CIM.leaf.QTL.Model)




## Fit a multiple-QTL model 
K21leafQTLs<- makeqtl(maprobR21,chr=c("1B","2D","3A"),pos=c(81,2,114),
                      qtl.name=c("K21.1B.SIM.leaf.QTL","K21.2D.SIM.leaf.QTL","R21.3A.CIM.leaf.QTL"),
                      what="prob")
K21leafQTL.model<- fitqtl(maprobR21,pheno.col=4,qtl=K21leafQTLs,formula=y~Q1*Q2*Q3,method="hk")
summary(K21leafQTL.model)

#### No Significant epistasis 




#### RALEIGH 21 HEADING DATE (JULIAN) ####

#column 6

#determine LOD threshold
out.permR21HDJulian<-scanone(maprobR21,pheno.col=6,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permR21HDJulian, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.50
#5%  3.16

#save LOD result
save(out.permR21HDJulian, file = "D:/HilliardxGA/QTL-RIL-results/LODperms/R21-HG-HDJulian-LODperm1000")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/HilliardxGA/QTL-RIL-results/LODperms/R21-HG-HDJulian-LODperm1000")

load("D:/HilliardxGA/QTL-RIL-results/R21-HD-SIM-results")
load("D:/HilliardxGA/QTL-RIL-results/R21-HD-CIM-results")

# SIM
SIMR21HD<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMR21HD<-scanone(SIMR21HD,pheno.col=6,model="normal",method="em")
summary(out.SIMR21HD,alpha=0.05,perms=out.permR21HDJulian,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr  pos ci.low ci.high lod  pval
#S5A_581137922  5A 79.4     15     190 3.3 0.033


# CIM
#window size 6; 4 covariate markers
out.CIMR21HD<-cim(SIMR21HD,pheno.col=6, n.marcovar=4,method="em",window=6)
summary(out.CIMR21HD,alpha=0.05,perms=out.permR21HDJulian,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr  pos ci.low ci.high  lod  pval
#S1D_425671465   1D 56.4   54.2      60 3.90 0.011
#S5A_584035850   5A 81.2   76.0      83 6.78 0.000
#c7B.1.loc17   7B.1 17.0   16.0      23 4.76 0.001

#save QTL results

save(out.SIMR21HD, file = "D:/HilliardxGA/QTL-RIL-results/R21-HD-SIM-results")
save(out.CIMR21HD, file = "D:/HilliardxGA/QTL-RIL-results/R21-HD-CIM-results")
write.table(out.SIMR21HD,file="D:/HilliardxGA/QTL-RIL-results/R21-HD_SIM.csv", sep = ",")
write.table(out.CIMR21HD,file="D:/HilliardxGA/QTL-RIL-results/R21-HD_SIM.csv", sep = ",")


#### MODEL PVE AND EPISTASIS ####

#First, use 'makeqtl' to declare the positions of putative QTL
#then, use fitqtl to model the effects of these QTL; for brevity, use the LSmeans to obtain a multienvironmental model


# a) SIM - Fit single QTL model


R21HD.5A.1.QTL<-makeqtl(maprobR21,chr="5A",pos=92.74202,qtl.name="S5A_536660424",what="prob")
R21HD.5A.1.QTL.Model<-fitqtl(maprobR21,pheno.col=6,qtl=R21HD.5A.1.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R21HD.5A.1.QTL.Model)

R20HD.5A.2.QTL<-makeqtl(maprobR21,chr="5A",pos=193,qtl.name="S5A_581137922",what="prob")
R20HD.5A.2.QTL.Model<-fitqtl(maprobR21,pheno.col=6,qtl=R20HD.5A.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.5A.2.QTL.Model)

## Fit a multiple-QTL model 
R21HDQTLs<- makeqtl(maprobR21,chr=c("5A","5A"),pos=c(92.74202,193),
                    qtl.name=c("R21HD.5A.1.QTL","R20HD.5A.2.QTL"),
                    what="prob")
R21HDMQTL.model<- fitqtl(maprobR21,pheno.col=6,qtl=R21HDQTLs,formula=y~Q1*Q2,method="hk")
summary(R21HDMQTL.model)


# b) out.CIM
## Fit single QTL model
R20HD.2B.QTL<-makeqtl(maprobR21,chr="2B",pos=88,qtl.name="c2B.loc88",what="prob")
R20HD.2B.QTL.Model<-fitqtl(maprobR21,pheno.col=10,qtl=R20HD.2B.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.2B.QTL.Model)

R20HD.5A.QTL<-makeqtl(maprobR21,chr="5A",pos=237,qtl.name="S5A_570227648",what="prob")
R20HD.5A.QTL.Model<-fitqtl(maprobR21,pheno.col=10,qtl=R20HD.5A.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(R20HD.5A.QTL.Model)

CIMR20HD.5B.2.QTL<-makeqtl(maprobR21,chr="5B.2",pos=4,qtl.name="c5B.2.loc4",what="prob")
CIMR20HD.5B.2.QTL.Model<-fitqtl(maprobR21,pheno.col=10,qtl=CIMR20HD.5B.2.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR20HD.5B.2.QTL.Model)

CIMR20HD.7D.QTL<-makeqtl(maprobR21,chr="7D",pos=16,qtl.name="c7D.loc16",what="prob")
CIMR20HD.7D.QTL.Model<-fitqtl(maprobR21,pheno.col=10,qtl=CIMR20HD.7D.QTL,formula=y~Q1,method="hk",get.ests=TRUE,dropone=TRUE)
summary(CIMR20HD.7D.QTL.Model)


## Fit a multiple-QTL model 
CIMR20HDQTLs<- makeqtl(LAmaprobR21,chr=c("2B","5A","5B.2","7D"),pos=c(88,237,4,16),
                    qtl.name=c("R20HD.2B.QTL","R20HD.5A.QTL","CIMR20HD.5B.2.QTL", "CIMR20HD.7D.QTL"),
                    what="prob")
CIMR20HDMQTL.model<- fitqtl(LAmaprobR21,pheno.col=10,qtl=CIMR20HDQTLs,formula=y~Q1*Q2*Q3*Q4,method="hk")
summary(CIMR20HDMQTL.model)


#### KINSTON 21 HEADING DATE (JULIAN) ####

#column 7

#determine LOD threshold
out.permK21HDJulian<-scanone(maprobR21,pheno.col=7,n.perm=1000, verbose=TRUE,method="em")#,model="normal")
summary(out.permK21HDJulian, alpha=c(0.20,0.05))
#LOD thresholds (1000 permutations)
#lod
#20% 2.49
#5%  3.18

#save LOD result
save(out.permK21HDJulian, file = "D:/HilliardxGA/QTL-RIL-results/LODperms/K21-HG-HDJulian-LODperm1000")

#after first time, quick load:
LAmaprobR21<- calc.genoprob(LAmap, step=1,map.function="kosambi",error.prob=0.001)#hidden Markov model
load("D:/HilliardxGA/QTL-RIL-results/LODperms/K21-HG-HDJulian-LODperm1000")

load("D:/HilliardxGA/QTL-RIL-results/K21-HD-SIM-results")
load("D:/HilliardxGA/QTL-RIL-results/K21-HD-CIM-results")

# SIM
SIMLAK20HD<-calc.genoprob(map, step=1,map.function="kosambi",error.prob=0.001)
out.SIMLAK20HD<-scanone(SIMLAK20HD,pheno.col=7,model="normal",method="em")
summary(out.SIMLAK20HD,alpha=0.05,perms=out.permK21HDJulian,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr  pos ci.low ci.high  lod  pval
#S5A_581137922   5A 79.4     70     189 4.06 0.005
#c7B.1.loc15   7B.1 15.0      6      23 4.42 0.004


# CIM
#window size 6; 4 covariate markers
out.CIMLAK20HD<-cim(SIMLAK20HD,pheno.col=7, n.marcovar=4,method="em",window=6)
summary(out.CIMLAK20HD,alpha=0.05,perms=out.permK21HDJulian,format="tabByCol",pvalues=TRUE)%T>% {print(summary(.))}
#    Length Class      Mode
#lod 6      data.frame list
#chr  pos ci.low ci.high  lod  pval
#S5A_581137922   5A 79.4     76      83 8.29 0.000
#c7B.1.loc17   7B.1 17.0     16      23 3.98 0.006


#save QTL results
save(out.SIMLAK20HD, file = "D:/HilliardxGA/QTL-RIL-results/K21-HD-SIM-results")
save(out.CIMLAK20HD, file = "D:/HilliardxGA/QTL-RIL-results/K21-HD-CIM-results")
write.table(out.SIMLAK20HD,file="D:/HilliardxGA/QTL-RIL-results/K21-HD_SIM.csv", sep = ",")
write.table(out.CIMLAK20HD,file="D:/HilliardxGA/QTL-RIL-results/K21-HD_CIM.csv", sep = ",")



#### MODEL PVE AND EPISTASIS ####

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










#### PLOT SIGNIFICANT QTL RESULTS ####

library('LinkageMapView')
#load map from file
setwd("C:/Users/Daniela/Documents/Wax/HGfinal/")
map<-read.cross(format="csv",file="UX1992-HilliardxGA-RIL-Miller-3Alinkage-map.csv",estimate.map=FALSE,
                na.strings=c(".","-","NA"),genotypes=c("A","B"),crosstype="riself")%>%jittermap(.)%>%calc.genoprob(map.function="kosambi")
summary(map)
setwd("C:/Users/Daniela/Documents/Wax/HGfinal/QTL-results/")

#DENSITY PLOT
lmv.linkage.plot(map, "HG_density.pdf", denmap = TRUE)

#MAKE QTL RESULTS DATAFRAME

# make a df to pass qtl info
qtldf <- data.frame(
  chr = character(),
  qtl = character(),
  so = numeric(),
  si = numeric(),
  ei = numeric(),
  eo = numeric(),
  col = character(),
  stringsAsFactors = FALSE
)
#ADD SIGNIFICANT QTL INDIVIDUALLY
#FOLLOW BELOW CONVENTION
# so = confidence interval low cM
# eo = confidence interval high cM
# si = same as so
# ei = same as eo
# FORESTGREEN color for LEAF 
# DARKORANGE color for SPIKE

#leaf QTL
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "1A",
                 qtl = "QFlg.ncb-1A",
                 so = 80.0,
                 si = 80.0,
                 ei = 92.5,
                 eo = 92.5,
                 col="forestgreen"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "2A",
                 qtl = "QFlg.ncb-2A",
                 so = 0.0,
                 si = 0.0,
                 ei = 1.5,
                 eo = 1.5,
                 col="forestgreen"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "2D",
                 qtl = "QFlg.ncb-2D",
                 so = 0.0,
                 si = 0.0,
                 ei = 1.6,
                 eo = 1.6,
                 col="forestgreen"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "3A",
                 qtl = "QFlg.ncb-3A",
                 so = 119.0,
                 si = 119.0,
                 ei = 122,
                 eo = 122,
                 col="forestgreen"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "4B",
                 qtl = "QFlg.ncb-4B",
                 so = 17,
                 si = 17,
                 ei = 20,
                 eo = 20,
                 col="forestgreen"
               ))

#spike QTL
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "1A",
                 qtl = "QSpg.ncb-1A.1",
                 so = 0,
                 si = 0,
                 ei = 1.6,
                 eo = 1.6,
                 col="darkorange"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "3A",
                 qtl = "QSpg.ncb-3A",
                 so = 122,
                 si = 122,
                 ei = 125,
                 eo = 125,
                 col="darkorange"
               ))
qtldf <- rbind(qtldf,
               data.frame(
                 chr = "5B",
                 qtl = "QSpg.ncb-5B",
                 so = 36,
                 si = 36,
                 ei = 39,
                 eo = 39,
                 col="darkorange"
               ))

#best version wide
lmv.linkage.plot(map, "HG_QTL_plot.pdf", mapthese = c("1A","2A","2D","3A","4B","5B"),
                 lgw = 0.2, pdf.height = 8, pdf.width = 12,
                 lgperrow = 3, lcol = "darkgray", rcol = "black", rfont=3, lfont=2,labdist = 0.1,
                 qtldf=qtldf, lg.col = "beige", rcex = 1, ruler = TRUE,
                 showonly=c("S1A_6370170","S1A_535882457","S2A_1872142","S2A_771789729","S2D_22461473",
                            "S3A_1287675","S3A_737432147","S4B_132328665","S4B_662157150","S5B_14987934","S5B_641803649",
                            "S2A_1872142","S2D_17164152","S4B_559855530","S1A_6370170","S3A_621010472","S5B_78423903",
                            "S1A_520868230","S1A_531680970","S3A_613343859","S3A_616867989","S1A_8516115",
                            "S1A_512491739","S2A_13121298","S2D_17164152","S2D_18184603",
                            "S4B_541252759","S4B_565312655","S5B_87491939","S5B_81171755"
                 ))


#3A chrom only

#add marker formatting
# make a list to pass label options
flist <- list()
#bold shared markers
locus <- c("S3A_567971108","S3A_609909640","S3A_613343859","S3A_621010472","S3A_627915876")
font  <- c(4)   #bold italic
flist[[1]] <- list(locus = locus, font = font)
#font  <- c(4)   #bold italic c(2) bold
locus <- c("S3A_613343859", "S3A_609909640")
font  <- c(4)   #bold italic
col <- c("red")
flist[[2]] <- list(locus = locus, font = font, col = col)
locus <- c("S3A_616867989")
col <- c("red")
flist[[3]] <- list(locus = locus, col = col)

#3A plot
#best version
lmv.linkage.plot(map, "HG_QTL_3A.pdf", mapthese = c("3A"),
                 lgperrow = 8, lcol = "darkgray", rcol = "black", rfont=3, lfont=1,
                 lg.col = "beige", rcex = 1, ruler = FALSE, pdf.height = 5,
                 markerformatlist = flist, labdist = 0.1,
                 showonly=c("S3A_1287675","S3A_737432147","S3A_613343859","S3A_616867989",
                            "S3A_609909640","S3A_579373891","S3A_617340619",
                            "S3A_621010472","S3A_627915876","S3A_567971108"
                 ))
