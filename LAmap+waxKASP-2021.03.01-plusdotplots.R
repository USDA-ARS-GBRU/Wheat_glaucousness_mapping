##### LA95135 x ASG2000 RIL population #####
#linkage map construction
#original map & parameters from Dr. Luis Rivera-Burgos
#i added KASP markers around 2A& 3A wax QTL

library("qtl")
library("magrittr")
library("ASMap")


#format my wax KASP data to match map
## encode alleles to match previous markers

#A allele is from LA95135 parent
#B allele is from AGS2000 parent



######### FIRST SECTION TO MERGE KASP MARKERS ONLY; NOT BUILDING MAP YET ######

## read in my new KASP markers as csv to fix formatting first
setwd("D:/wax/LA-KASP-fullpop")
KASPw <- read.csv(file="waxKASP_NCB181289A-1292A.csv")
head(KASPw)

#create correct name column to match map
KASPw$genotype <- paste(KASPw$X, KASPw$X.1, sep="")
head(KASPw)
#write out new df
write.csv(KASPw, "waxKASP_NCB181289A-1292A.csv")


#merge KASP data files - my wax KASP markers + KASP from Luis (previous)

## read in KASP from Luis
setwd("D:/wax/LA-linkage-map/raw-from-Luis")
KASPL <- read.csv(file="KASP.6.19.2020.csv")
head(KASPL)
ncol(KASPL)

## read in my KASP
setwd("D:/wax/LA-KASP-fullpop")
KASPw <- read.csv(file="waxKASP_NCB181289A-1292A_mapformat.csv")
head(KASPw)
ncol(KASPw)

#merge KASP files
KASPall <- merge(x=KASPL, y=KASPw, by = "genotype", all.x = TRUE)
head(KASPall)
ncol(KASPall)

#write out all KASP markers to one file
setwd("D:/wax/LA-linkage-map")
write.csv(KASPall,file="KASP-all-2021.02.22.csv")


##### Read in data #####

setwd("D:/wax/LA-linkage-map")

LApop_geno<-read.cross(format="csv",file="LApopGeno26LG_6.19.2020.csv",estimate.map=FALSE,
                       na.strings=c("-","NA"),genotypes=c("A","B"),
                       crosstype="riself")%T>%{print(summary(.))}

#A allele is from LA95135 parent
#B allele is from AGS2000 parent

summary(LApop_geno)


## read in KASP markers

KASP<-read.cross(format="csv",file="KASP-all-2021.02.22.csv",estimate.map=FALSE,
                 na.strings=c("NA"),genotypes=c("A","B"),
                 crosstype="riself")%T>% {print(summary(.))}
summary(KASP)







##### BELOW DIRECTLY FROM LUIS CODE #####

#-Pre-Construction
#***************************************************************************************************************************************************
# Checking GENOTYPES with large ammount of missing marker data
#plotMissing(LApop_geno) #check for darkerst lines meaning 
sg<-statGen(LApop_geno, bychr = FALSE, stat.type =c("miss"), id="genotype")%T>% {print(.)}
LApop_geno1<-subset(LApop_geno, ind = sg$miss <6000) %T>% {print(summary(.))}#6000 is 50% of marker missing data per genetoype
#-From a map construction point of view, highly related individuals can enhance segregation distortion of markers
#gc<-genClones(LApop_geno1,tol=0.95,id="genotype") %T>% {print(.$cgd)}
#*** No problems at this level. We keep all Genotypes adn the same working file

#-Check MARKERS for segregation distortion
#The level of segregation distortion, the allelic proportions and the missing value proportion
#across the genome can be graphically represented using the marker profiling function
#profileMark()
#options(max.print=13000)
#seg_dis<-profileMark(LApop_geno1,stat.type=c("seg.dist","prop","miss"),crit.val="bonf",layout=c(1,4),type="p",cex=0.5)# %T>% {print(.)}
#setwd("C:/Users/lariver2/Documents/LA95135xAGS2000/LApop_Map_&_QTLanalysis/ASMap/LAMAP_plus_KASP") #set working directory
#write.csv(seg_dis,file="LApop_SNPs_Seg_Dis.6.19.20.csv")

#mm<-statMark(LApop_geno1,stat.type="marker")$marker$AA %T>%{print(summary(.))}
#pm<-profileMark(LApop_geno1,stat.type=c("seg.dist","dxo","erf","lod"),id ="genotype",layout=c(1, 4),type="l")%T>% {print(summary(.))}

#-Droping marker with high segregation distortion and missing data
LApop_geno2<-pullCross(LApop_geno1,type="missing",pars=list(miss.thresh =0.1))
LApop_geno2$missing$table
LApop_geno3<-pullCross(LApop_geno2,type="seg.distortion",pars=list(seg.thresh=0.01))#%T>% {print(summary(.))}## The seg distorted markers were removed # 112377 SNPs
LApop_geno3$seg.distortion$table
LApop_geno4<-pullCross(LApop_geno3,type="co.located") #delete duplicated SNPs patterns
LApop_geno4$co.located$table
#Adding (ajouter) KASP markers
#KASP<-read.cross(format="csv",file="KASP.6.19.2020.csv",estimate.map=FALSE,
#                 na.strings=c("NA"),genotypes=c("A","B"),
#                 crosstype="riself")%T>% {print(summary(.))}
#KASP$geno
LApop_geno5<-combineMap(LApop_geno4,KASP,id="genotype",keep.all=TRUE,
                        merge.by = "genotype")
jittermap(LApop_geno5)
summary(LApop_geno5)
class(LApop_geno5)
#-Build the map
LAgeneticMap<- mstmap(LApop_geno5,pop.type="ARIL",dist.fun="kosambi",id="genotype",
                      objective.fun="COUNT",noMap.dist=15,noMap.size=0,bychr=TRUE,
                      miss.thresh = 1,mvest.bc = FALSE,detectBadData = FALSE, 
                      as.cross=TRUE, return.imputed=FALSE,trace=TRUE,anchor=TRUE,
                      p.value=1e-12)

nmar(LAgeneticMap)
chrlen(LAgeneticMap)
summaryMap<-summary.map(LAgeneticMap)%T>%{print(.)}
#NOTE: if chromosomes of D subgenome have very few markers, join them in sub-groups so they dont get dropped in the next step
#Example: MMCross<-mergeCross(MMCross,merge=list("4B"=c("4B.1","4B.2")))
#NOTE: remove all linkage groups with 10 or fewer markers
chrnamesDROP<-rownames(subset(summaryMap,summaryMap$n.mar <10))
markernamesDROP<-markernames(LAgeneticMap,chr=chrnamesDROP)
LAgeneticMap<-drop.markers(LAgeneticMap,markers=markernamesDROP)
summary.map(LAgeneticMap)
##CHECK if the marker pos (bp) is sorted from lower to higher within each LG
#markernames(LAgeneticMap,chr="1A.10")
#-Merge subgroups to get to 26 linkage groups
#Note.- You can make 21, 25 or 26 LG by modifiying "mergeCross" command
#Method 2.- Join only cluster with maximun number of SNP-markers
Merged_LAgeneticMap<-mergeCross(LAgeneticMap,merge=list
                                ("1D"=c("1D.1","1D.4"),
                                  "2D"=c("2D.1","2D.2","2D.4"),
                                  "3A"=c("3A.1","3A.3"),
                                  "5A"=c("5A.1","5A.4"),
                                  "6B"=c("6B.1","6B.2"),
                                  "7B"=c("7B.2","7B.6"),"7D"=c("7D.4","7D.1")))
summary.map(Merged_LAgeneticMap)
#heatMap(Merged_LAgeneticMap)
#-Reconstructing the Genetic Map
Final_LAgeneticMap<-mstmap(Merged_LAgeneticMap,pop.type="ARIL",dist.fun="kosambi",
                           id="genotype",objective.fun="COUNT",noMap.dist=15,noMap.size=0,
                           bychr=TRUE,miss.thresh=1,mvest.bc=FALSE,detectBadData=FALSE, 
                           as.cross=TRUE, return.imputed=FALSE,trace=TRUE,anchor=TRUE,
                           p.value=1)
nmar(Final_LAgeneticMap)
#updated output:
#1A.1 1B.1 1B.4   1D 2A.2 2B.1   2D   3A 3B.1 3D.4 4A.1 4A.4 4B.1 4D.1   5A 5B.2 5B.4 5D.5   6A   6B 6D.1 6D.6 7A.1 7A.2   7B   7D 
#275   65   21  171   89  415  104  240  404   18   49   48  172   22  209  107   69   23  137  117   20   41  130  172  260   30 
jittermap(Final_LAgeneticMap)
chrlen(Final_LAgeneticMap)
#updated output:
#     1A.1      1B.1      1B.4        1D      2A.2      2B.1        2D        3A      3B.1      3D.4      4A.1      4A.4      4B.1      4D.1        5A      5B.2 
#279.76194  51.47373  40.30501 197.26031 105.40112 302.33491 242.32304 282.60334 358.74908  48.82355 129.54588  38.32316 174.90623  80.76779 341.85708 162.70678 
#5B.4      5D.5        6A        6B      6D.1      6D.6      7A.1      7A.2        7B        7D 
#40.08642  51.67796 168.61474 187.74299  46.17836  73.69567 127.73639 179.73496 245.03234 150.71106
names(Final_LAgeneticMap$geno)<-c("1A","1B.1","1B.2","1D","2A","2B","2D","3A","3B","3D",
                                  "4A.1","4A.2","4B","4D","5A","5B.2","5B.1","5D","6A",
                                  "6B","6D.1","6D.2","7A.1","7A.2","7B","7D")%T>%{print(.)}
summary.map(Final_LAgeneticMap)
plot.map(Final_LAgeneticMap)
## NOTE to myselef: chromosomes 5B.2 and 5B.1 are flipped as LG. The markers within LG
# are in the right order. I have to flip it in excel and make figures and QTL analysis.
#heatMap(Final_LAgeneticMap, main=NULL)
class(Final_LAgeneticMap)

#********** Save_the_Final_Map_as_CrossFile_for the QTL ANALYSIS *****************
setwd("D:/wax/Field2020/QTL-mapping_waxKASP")
write.cross(Final_LAgeneticMap,"csv",filestem="LAGenMap_26LG+waxKASP_2021.03.01",
            c("1A","1B.1","1B.2","1D","2A","2B","2D","3A","3B","3D",
              "4A.1","4A.2","4B","4D","5A","5B.2","5B.1","5D","6A",
              "6B","6D.1","6D.2","7A.1","7A.2","7B","7D"))

#LAsummary <- summary.map(Final_LAgeneticMap)
#LAsummary
#write.csv(LAsummary, "LAmap+waxKASP.csv")
#can't export this kinda object


#### THIS IS THE CURRENT FINAL MAP VERSION - FOR USE AND SHARE 


# FINAL CHECK DOT PLOT AFTER REMOVING OUTLIERS 
dotplot<-read.csv(file = "D:/wax/Field2020/QTL-mapping_waxKASP/dotplot_LAGenMap_26LG+waxKASP_2021.03.01.csv", header = TRUE)
head(dotplot,n=5)
tail(dotplot,n=5)

library(ggplot2)
library(grid)
ggplot(dotplot,aes(cM,bp))+geom_point(size=1,colour="black")+facet_wrap(~Chr, nrow=7)
#saved plot to file: UX1992-HilliardxGA-Miller-linkagemap-dotplots.jpg and pdf