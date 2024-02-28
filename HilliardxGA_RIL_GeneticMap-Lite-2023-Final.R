#-BUILDING A GENETIC MAP FOR DH POPULATION- MD233xSS8641
#install.packages("onemap")


#first, vcf to R/QTL format:

  ## Only need to run lines below once to install
install.packages("devtools")
library("devtools")
#install_github("youruser/yourrepo")
install_github("etnite/bwardr",force = TRUE)
  
# Open libraries in R
library(devtools)
library(gaston)
library(bwardr); packageVersion("bwardr")

#### INITIAL VCF FILE CONVERSION ####


setwd("D:/HilliardxGA/")

bed <- read.vcf(file="SunRILs_2021_UX1992_RILs_filt.vcf.gz", convert.chr = FALSE)

abh <- format_qtlmap_geno(bed = bed,
                          par_a = "HILLIARD",
                          par_b = "GA06493-13LE6",
                          rm_het = TRUE, #Note: when FALSE show SNPs parents H calls.
                          #When true removes SNPs parents H calls
                          rm_miss = TRUE,
                          include_pars = TRUE,
                          out_fmt = "rqtl")


write.csv(abh$abh,"HilliardxGA06493-13LE6_filt_vcf_rqtl_formatted.csv", row.names = FALSE)
#write.csv(abh,"HilliardxGA06493-13LE6_filt_vcf_rqtl_formatted.csv", row.names = FALSE)


############ ALLELE DEFINITIONS ##############

###
# par_a = "HILLIARD" ---> PARENT ALLELE "A"
# par_b = "GA06493-13LE6" ---> PARENT ALLELE "B"





########### BUILD LINKAGE MAP ############


install.packages("qtl")
install.packages("magrittr")
install.packages("ASMap")
#options(max.print=999999)
memory.limit(size = NA) #64 GB RAM max [USDA PC]
library("qtl"); packageVersion("qtl")
#version '1.48.1'
library("magrittr"); packageVersion("magrittr")
#version '2.0.1'
library("ASMap"); packageVersion("ASMap")
#version '1.0.4'
DH_geno<-read.cross(format="csv",dir="D:/HilliardxGA/",
                    file="HilliardxGA06493-13LE6_filt_vcf_rqtl_formatted_transpose_noparents.csv",estimate.map=FALSE,na.strings=c("-","NA"),genotypes=c("A","B"),
                    crosstype="riself") %T>% {print(summary(.))}

#FILTERING MISSING DATA & SEGREGATION DISTORTION 
sg<-statGen(DH_geno, bychr = FALSE, stat.type =c("miss"), id="genotype")%T>% {print(.)}
DH_geno1<-subset(DH_geno, ind = sg$miss <(10347/2)) %T>% {print(summary(.))} #filter individuals with >50% missing data given 10347 starting markers
DH_geno2<-pullCross(DH_geno1,type="missing",pars=list(miss.thresh =0.10))%T>% {print(summary(.))} #2534 markers remaining
DH_geno3<-pullCross(DH_geno2,type="seg.distortion",pars=list(seg.thresh=0.01))%T>% {print(summary(.))}# 2368 markers remaining
DH_geno4<-pullCross(DH_geno3,type="co.located") %T>% {print(summary(.))} #delete duplicated SNPs patterns; 1538 markers remaining
summary(DH_geno4)

#BUILD THE GENETIC MAP
DHMap<- mstmap(DH_geno4,pop.type="riself",dist.fun="kosambi",id="genotype",objective.fun="COUNT", 
                      noMap.dist=15,noMap.size = 0,bychr=TRUE,miss.thresh = 1,mvest.bc = FALSE,detectBadData = FALSE, 
                      as.cross=TRUE, return.imputed=FALSE,trace=TRUE,anchor=TRUE,p.value=1e-16)
summaryMap<-summary.map(DHMap)%T>%{print(.)}



#FILTER LINKAGE GROUPS WITH LOW MARKER COVERAGE

#if chromosomes of D subgenome have very few markers, join them in sub-groups so they dont get dropped in the next step
#DHMap<-mergeCross(DHMap,merge=list("2D"=c("2D.2","2D.3"),
#                                   "6D"=c("6D.1","6D.5"),"7D"=c("7D.1","7D.2","7D.3","7D.4","7D.5","7D.6")))
#summaryMap<-summary.map(DHMap)%T>%{print(.)}


#remove all linkage groups with 10 or fewer markers
chrnamesDROP<-rownames(subset(summaryMap,summaryMap$n.mar <10)) #keep D genome LGs with only 8 markers... 4D only has 3 though :'(
markernamesDROP<-markernames(DHMap,chr=chrnamesDROP)
DHMap<-drop.markers(DHMap,markers=markernamesDROP)
summary.map(DHMap)
##CHECK if the marker pos (bp) is sorted from lower to higher within each LG
View(DHMap)

#FLIP LINKAGE GROUPS NOT IN ASCENDING ORDER
DHMap_3A_1_flip <- flip.order(DHMap, "3A.1")
#View(DHMap_3A_1_flip)
DHMap_ordered <- flip.order(DHMap_3A_1_flip, "7B.1")
#View(DHMap_ordered)
plot.map(DHMap_ordered)
summary.map(DHMap_ordered)

#MERGE SUBGROUPS TO GET 26 LINKAGE GROUPS (IF POSSIBLE--AVOID LGS W >30 cM GAPS)
Merged_DHMap<-mergeCross(DHMap_ordered,merge=list
                         ("1B" = c("1B.1","1B.3"),
                          "2A" = c("2A.1","2A.2"),
                          "2B" = c("2B.1","2B.2"),
                          "3A" = c("3A.1","3A.2","3A.3"),
                          "3B" = c("3B.1","3B.6"),
                          "4A" = c("4A.1","4A.3"),
                          "5A" = c("5A.3","5A.5")
                           ))
summary.map(Merged_DHMap)

#drop 6A.2
markersDROP<-markernames(Merged_DHMap,chr="6A.2")
Finale_DHMap<-drop.markers(Merged_DHMap,markers=markersDROP)
summary.map(Finale_DHMap)


#RECONSTRUCTING THE GENETIC MAP AFTER MERGING LINKAGE GROUPS
Final_DHMap<-mstmap(Finale_DHMap,pop.type="riself",dist.fun="kosambi",id="genotype",objective.fun="COUNT", 
                    noMap.dist=15,noMap.size = 0,bychr=TRUE,miss.thresh = 1,mvest.bc = FALSE,detectBadData = FALSE, 
                    as.cross=TRUE, return.imputed=FALSE,trace=TRUE,anchor=TRUE,p.value=1)
names(Final_DHMap$geno)<-c("1A","1B","1D","2A","2B","2D","3A","3B",
                           "4A","4B","5A","5B","6A","6B","6D",
                           "7A","7B.1","7B.2","7D")%T>% {print(.)}
summary.map(Final_DHMap)
plot.map(Final_DHMap)
View(Final_DHMap)
#CHECK if bps within linkage groups are in ascending order



#SAVE MAP AS CROSS FILE FOR QTL ANALYSIS
setwd("D:/HilliardxGA")
write.cross(Final_DHMap,"csv",filestem="HilliardxGA-RIL-linkage-map",
            c("1A","1B","1D","2A","2B","2D","3A","3B",
              "4A","4B","5A","5B","6A","6B","6D",
              "7A","7B.1","7B.2","7D"))




#CREATE DOTPLOT TO CHECK MARKER ORDERING -- FOR GBS MARKERS ONLY (HAVE BP DATA)
dotplot<-read.csv(file = "D:/HilliardxGA/dotplot-UX1992-linkage-map.csv", header = TRUE)
head(dotplot,n=5)
tail(dotplot,n=5)

library(ggplot2)
library(grid)
ggplot(dotplot,aes(cM,bp))+geom_point(size=1,colour="black")+facet_wrap(~Chr, nrow=7)
#saved plot to file: UX1992-HilliardxGA-Miller-linkagemap-dotplots-semifinal.jpg and pdf


#REMOVE OUTLIER MARKERS IDENTIFIED FROM ABOVE DOTPLOT 
outliers <- c("S1A_527395289", "S1A_208718176", "S1B_393838944")
Final_Map_UX1992<-drop.markers(Final_DHMap,markers=outliers)
summary.map(Final_Map_UX1992)


#SAVE MAP AS CROSS FILE FOR QTL ANALYSIS
setwd("D:/HilliardxGA")
write.cross(Final_Map_UX1992,"csv",filestem="UX1992-HilliardxGA-RIL-Miller-linkage-map",
            c("1A","1B","1D","2A","2B","2D","3A","3B",
              "4A","4B","5A","5B","6A","6B","6D",
              "7A","7B.1","7B.2","7D"))
#### THIS IS THE CURRENT FINAL MAP VERSION - FOR USE AND SHARE !



# SUMMARY AND PLOT OF FINAL LINKAGE MAP 
summary.map(Final_Map_UX1992)
plot.map(Final_Map_UX1992)



# FINAL CHECK DOT PLOT AFTER REMOVING OUTLIERS 
dotplot<-read.csv(file = "D:/HilliardxGA/dotplot-UX1992-linkage-map-clean.csv", header = TRUE)
head(dotplot,n=5)
tail(dotplot,n=5)

library(ggplot2)
library(grid)
ggplot(dotplot,aes(cM,bp))+geom_point(size=1,colour="black")+facet_wrap(~Chr, nrow=7)
#saved plot to file: UX1992-HilliardxGA-Miller-linkagemap-dotplots.jpg and pdf
