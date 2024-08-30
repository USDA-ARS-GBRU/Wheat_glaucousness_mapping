## Entry means calculation for wax ratings
# HilliardxGA-13LE6 

##### load libraries #####

#install emmeans package
#library(devtools)
#packageVersion('devtools') #must be v2.0 or higher
#remotes::install_github("rvlenth/emmeans")

library("lme4")
library("lmerTest")
library("emmeans")

#### IMPORT data files #############

### RALEIGH all wax trts R21 LA pop
R21<-read.table("D:/wax/Field2021/final-export-2021fieldbook/R21_T20-T25_UX1992_final_DMM_Rready.csv",sep=",",head=TRUE,
                colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
dR21 <- as.Date(R21$flowering, format="%m/%d/%Y") #raw formatted in month-day-year -- default format: year-month-day
JulianDay21 <- julian(dR21, origin = as.Date("2021-01-01")) #set origin to first of the year from data
R21$JulianHD <- JulianDay21
R21$Wax.Spike <- R21$Wax.Spike.5.02.21
R21$Wax.Leaf <- R21$Wax.Leaf.Bottom.4.26.21
R21$Location <- "Raleigh, NC"
R21 <- subset (R21, select = c(-Wax.Spike.5.02.21, -Wax.Leaf.Bottom.4.26.21, -flowering))

### KINSTON all wax trts K21 HG pop
K21<-read.table("D:/wax/Field2021/final-export-2021fieldbook/K21-WaxExp_T6-T10_final_DMM_HGpop_Rready.csv",sep=",",head=TRUE,
                colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
#set HD col to "Date" class
dK21 <- as.Date(K21$flowering, format="%m/%d/%Y") #raw formatted in month-day-year
#convert "Date" class to "numeric" Julian date 
JulianDayK21 <- julian(dK21, origin = as.Date("2021-01-01")) #set origin to first of the year from data
#add Julian date col to dataframe
K21$JulianHD <- JulianDayK21
#rename columns
K21$Wax.Spike <- K21$Wax.Spike.4.22.21
K21$Wax.Leaf <- K21$Wax.Leaf.Bottom.4.22.21
#add location column
K21$Location <- "Kinston, NC"
#remove old names
K21 <- subset (K21, select = c(-Wax.Spike.4.22.21, -Wax.Leaf.Bottom.4.22.21, -flowering))

# concatenate RALEIGH and KINSTON datasets for combined environment means
total <- rbind(R21, K21)
#make sure rep is a factor variable
total$Rep <- as.factor(total$Rep)


#### spike means combined env model ####

#combined environment model
require(lme4)
mm_spike7 <- lmer(Wax.Spike ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, total)
summary(mm_spike7) 
isSingular(mm_spike7) #NOT SINGULAR
anova(mm_spike7)

#calculate means
emspike <- emmeans(mm_spike7,"Entry")
print(emspike)

#write output
setwd("~/Wax/HGfinal/means/")
write.csv(emspike,file="spike-emmeans-across-envs.csv")


#### leaf means combined env model ####

# combined env model
require(lme4)
mm_leafg <- lmer(Wax.Leaf ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, total)
isSingular(mm_leafg)
anova(mm_leafg)

#calculate means 
emleaf <- emmeans(mm_leafg, "Entry")
print(emleaf)

#write output
write.csv(emleaf,file="leaf-emmeans-across-envs.csv")


#### R21 WAX LEAF #### 

#ANOVA
ANOVA_leafB<-lmer(Wax.Leaf~Entry+(1|Rep),REML=TRUE,R21)
anova(ANOVA_leafB)
summary(ANOVA_leafB)

#lsmeans
leafBmeansbygeno<-lsmeans(ANOVA_leafB,"Entry")
print(leafBmeansbygeno)

#write output
write.csv(leafBmeansbygeno,file="R21-HG-wax-leaf-bottom-lsmeans.csv")


#### R21 WAX SPIKE #### 

#ANOVA
ANOVA_spike<-lmer(Wax.Spike~Entry+(1|Rep),REML=TRUE,R21)
anova(ANOVA_spike)
summary(ANOVA_spike)

#lsmeans
spikemeansbygeno<-lsmeans(ANOVA_spike,"Entry")
print(spikemeansbygeno)

#write output
write.csv(spikemeansbygeno,file="R21-HG-wax-spike-lsmeans.csv")


#### K21 WAX LEAF #### 

#ANOVA
ANOVA_leafBK<-lmer(Wax.Leaf~Entry+(1|Rep),REML=TRUE,K21)
anova(ANOVA_leafBK)
summary(ANOVA_leafBK)

#lsmeans
leafBmeansbygenoK<-lsmeans(ANOVA_leafBK,"Entry")
print(leafBmeansbygenoK)

#write output
write.csv(leafBmeansbygenoK,file="K21-HG-wax-leaf-bottom-lsmeans.csv")


#### K21 WAX SPIKE #### 

#ANOVA
ANOVA_spikeK<-lmer(Wax.Spike~Entry+(1|Rep),REML=TRUE,K21)
anova(ANOVA_spikeK)
summary(ANOVA_spikeK)

#lsmeans
spikemeansbygenoK<-lsmeans(ANOVA_spikeK,"Entry")
print(spikemeansbygenoK)

#write output
write.csv(spikemeansbygenoK,file="K21-HG-wax-spike-lsmeans.csv")


