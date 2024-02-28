## Lsmeans calculation for wax leaf types
#Raleigh 2020 LA population

library("lme4")
#library("lmerTest")
#library("lsmeans")
library("emmeans")

#####
### all wax trts R20 LA pop
R21<-read.table("D:/wax/Field2021/final-export-2021fieldbook/R21_T20-T25_UX1992_final_DMM_Rready.csv",sep=",",head=TRUE,
                     colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
head(R21, n=5)
tail(R21, n=5)
colnames(R21)

#Julian Date

#set HD col to "Date" class
dR21 <- as.Date(R21$flowering, format="%m/%d/%Y") #raw formatted in month-day-year
#default format: year-month-day

#convert "Date" class to "numeric" Julian date 
JulianDay21 <- julian(dR21, origin = as.Date("2021-01-01")) #set origin to first of the year from data
class(JulianDay21)


#add Julian date col to dataframe
R21$JulianHD <- JulianDay21
colnames(R21)

R21$Wax.Spike <- R21$Wax.Spike.5.02.21
R21$Wax.Leaf <- R21$Wax.Leaf.Bottom.4.26.21

R21 <- subset (R21, select = c(-Wax.Spike.5.02.21, -Wax.Leaf.Bottom.4.26.21, -flowering))

#R21$Year <- "2021"
R21$Location <- "Raleigh, NC"
colnames(R21)
head(R21)



### KINSTON 2021 IMPORT 


### all wax trts K21 HG pop
K21<-read.table("D:/wax/Field2021/final-export-2021fieldbook/K21-WaxExp_T6-T10_final_DMM_HGpop_Rready.csv",sep=",",head=TRUE,
                colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
head(K21, n=5)
tail(K21, n=5)
colnames(K21)

#set HD col to "Date" class
dK21 <- as.Date(K21$flowering, format="%m/%d/%Y") #raw formatted in month-day-year
#default format: year-month-day

#convert "Date" class to "numeric" Julian date 
JulianDayK21 <- julian(dK21, origin = as.Date("2021-01-01")) #set origin to first of the year from data
class(JulianDayK21)


#add Julian date col to dataframe
K21$JulianHD <- JulianDayK21
head(K21)
colnames(K21)

K21$Wax.Spike <- K21$Wax.Spike.4.22.21
K21$Wax.Leaf <- K21$Wax.Leaf.Bottom.4.22.21

K21 <- subset (K21, select = c(-Wax.Spike.4.22.21, -flowering))


#K21$Year <- "2021"
K21$Location <- "Kinston, NC"

colnames(K21)
head(K21)
head(R21)


### COMBINE BOTH LOCATIONS IN ONE DATASET

all <- merge(R21, K21, by = "Entry")
?merge
all <- inner_join(R21, K21, by = "Entry")
head(all)
tail(all)
colnames(all)



total <- rbind(R21, K21)
colnames(total)
head(total)
tail(total)


########## SIMPLIFIED IMPORT #############

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


# concatenate RALEIGH and KINSTON datasets
total <- rbind(R21, K21)
write.csv(total, "D:/Coursework/2022-spring/ST537-Multivariate-Longitudinal-Data-Analysis/FinalProject/HillardxGA-R21K-combined-pheno.csv")

#######



ANOVA_spike <- lmer(Wax.Spike ~ Entry + (1|Location:Rep), REML=TRUE, total) #so far so good

#full model
fullmodel <- lmer(Wax.Spike ~ Entry + (1|Location) + (1|Location:Rep) + (1|JulianHD), REML=TRUE, total)
summary(fullmodel)
#full model with genotype "Entry" random effect - use output to look at relative contributions to variance
fullmodelr <- lmer(Wax.Spike ~ (1|Entry) + (1|Location) + (1|Location:Rep) + (1|JulianHD), REML=TRUE, total)
summary(fullmodelr)

#also try full model with heading date as a FIXED covariate - compare to other full model by AIC/BIC?
fullmodelf <- lmer(Wax.Spike ~ Entry + (1|Location) + (1|Location:Rep) + JulianHD, REML=TRUE, total)
summary(fullmodelf)

# reduced model without heading date covariate
noHD <- lmer(Wax.Spike ~ Entry + (1|Location) + (1|Location:Rep), REML=TRUE, total)
summary(noHD)



## This uses lme() not lmer() ...

# compare reduced model to full model - LRT _ pg29 in notes
anova.lme(fullmodel,noHD)



#failed to converge
ANOVA_spike <- lmer(Wax.Spike ~ Entry + (Wax.Spike|Location) + (Wax.Spike| Location:Rep), REML=TRUE, total)

# Mixed model - random effects: location, location:rep, heading date
ANOVA_spike <- lmer(Wax.Spike ~ 1 + Entry + (1|Location) + (1|Location:Rep) + (1|JulianHD), REML=TRUE, total)
anova(ANOVA_spike)
summary(ANOVA_spike)

#lsmeans
spikelsmeansbygeno<-lsmeans(ANOVA_spike,"Entry")
print(spikelsmeansbygeno)

#lsmeans
#(LAwaxheadsHD.rg <- ref.grid(ANOVA_spike1))
#LAmeansgenoHD<-lsmeans(ANOVA_spike1,"Entry")
#print(LAmeansgenoHD)

#write output
setwd("D:/HilliardxGA/means/")
write.csv(spikelsmeansbygeno,file="spike-lsmeans-all2021locs-HD.csv")



### model WITHOUT heading date effect

# Mixed model - random effects: location, location:rep
ANOVA_spike <- lmer(Wax.Spike ~ 1 + Entry + (1|Location) + (1|Location:Rep), REML=TRUE, total)
anova(ANOVA_spike)
summary(ANOVA_spike)
#lsmeans
spikelsmeansbygeno<-lsmeans(ANOVA_spike,"Entry")
#print(spikelsmeansbygeno)
#save
setwd("D:/HilliardxGA/means/")
write.csv(spikelsmeansbygeno,file="spike-lsmeans-all2021locs.csv")



#### LEAF means - combined model ####

# Mixed model - random effects: location, location:rep
ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|Location) + (1|Location:Rep), REML=TRUE, total)
anova(ANOVA)
summary(ANOVA)
entry_means<-lsmeans(ANOVA,"Entry")
setwd("D:/HilliardxGA/means/")
write.csv(entry_means,file="leaf-bottom-lsmeans-all2021locs.csv")

# Mixed model - random effects: location, location:rep, heading date
ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|Location) + (1|Location:Rep) + (1|JulianHD), REML=TRUE, total) # is singular 
anova(ANOVA)
#summary(ANOVA)
#lsmeans
entry_means<-lsmeans(ANOVA,"Entry")
#print(spikelsmeansbygeno)
#save
setwd("D:/HilliardxGA/means/")
write.csv(entry_means,file="leaf-bottom-lsmeans-all2021locs-HD.csv")

ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|JulianHD), REML=TRUE, total)
ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|Location:JulianHD), REML=TRUE, total)
ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|Location) + (1|Location:Rep) + (1|Location:JulianHD), REML=TRUE, total) # is singular 
ANOVA <- lmer(Wax.Leaf ~ 1 + Entry + (1|Location) + (1|Location:Rep) + (1|Location:JulianHD), REML=TRUE, total) # is singular 
anova(ANOVA)



#### wAX SPIKE ####

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
#?lmer
ANOVA_spike<-lmer(Wax.Spike.5.02.21~Entry+(1|Rep),REML=TRUE,R21)
anova(ANOVA_spike)
#?anova
summary(ANOVA_spike)

#trying to add JulianHD as a covariate
#ANOVA_spike1<-lmer(Wax.Spike.4.26.21~Entry+(1|Rep)+(1|JulianHD),REML=TRUE,R21)
#?isSingular
#isSingular(ANOVA_spike1)
#isSingular(ANOVA_spike)
#rePCA(ANOVA_spike1)
#anova(ANOVA_spike1, type="III")
#summary(ANOVA_spike1)

#lsmeans
#(waxheadsR.rg <- ref.grid(ANOVA_spike))
spikelsmeansbygeno<-lsmeans(ANOVA_spike,"Entry")
#Warning: Cannot use mode = "kenward-roger" because *pbkrtest* package is not installed
print(spikelsmeansbygeno)

#lsmeans
#(LAwaxheadsHD.rg <- ref.grid(ANOVA_spike1))
#LAmeansgenoHD<-lsmeans(ANOVA_spike1,"Entry")
#print(LAmeansgenoHD)

#write output
setwd("D:/wax/Field2021/GHpop-lsmeans")
write.csv(spikelsmeansbygeno,file="R21-GH-wax-spike-lsmeans.csv")




#### WAX LEAF BOTTOM #### 

#ANOVA
ANOVA_leafB<-lmer(Wax.Leaf.Bottom.4.26.21~Entry+(1|Rep),REML=TRUE,R21)
anova(ANOVA_leafB)
summary(ANOVA_leafB)

#lsmeans
#(LAleafB.rg <- ref.grid(ANOVA_leafB))
leafBmeansbygeno<-lsmeans(ANOVA_leafB,"Entry")
print(leafBmeansbygeno)

#write output
write.csv(leafBmeansbygeno,file="R21-HG-wax-leaf-bottom-lsmeans.csv")





#### HEADING DATE - JULIAN DAYS ####

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
ANOVA_HD<-lmer(JulianHD~Entry+(1|Rep),REML=TRUE,R21)
anova(ANOVA_HD)
summary(ANOVA_HD)

#lsmeans
#(LAHDR.rg <- ref.grid(ANOVA_HD))
HDlsmeansbygeno<-lsmeans(ANOVA_HD,"Entry")
print(HDlsmeansbygeno)

#write output
write.csv(HDlsmeansbygeno,file="R21-HG-HDjulian-lsmeans.csv")





#### 2021 KINSTON DATA #####

### all wax trts K21 HG pop
K21<-read.table("D:/wax/Field2021/final-export-2021fieldbook/K21-WaxExp_T6-T10_final_DMM_HGpop_Rready.csv",sep=",",head=TRUE,
                colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
head(K21, n=5)
tail(K21, n=5)
colnames(K21)

#set HD col to "Date" class
dK21 <- as.Date(K21$flowering, format="%m/%d/%Y") #raw formatted in month-day-year
#default format: year-month-day

#convert "Date" class to "numeric" Julian date 
JulianDayK21 <- julian(dK21, origin = as.Date("2021-01-01")) #set origin to first of the year from data
class(JulianDayK21)


#add Julian date col to dataframe
K21$JulianHD <- JulianDayK21
head(K21)
colnames(K21)


#### K21 HEADING DATE ########

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
ANOVA_HD21K<-lmer(JulianHD~Entry+(1|Rep),REML=TRUE,K21)
anova(ANOVA_HD21K)
summary(ANOVA_HD21K)

#lsmeans
#(LAHDR.rg <- ref.grid(ANOVA_HD20))
HDlsmeans21K<-lsmeans(ANOVA_HD21K,"Entry")
print(HDlsmeans21K)

#write output
write.csv(HDlsmeans21K,file="K21-HG-HDjulian-lsmeans.csv")


## K21 wAX SPIKE 

#ANOVA
ANOVA_spike_K21<-lmer(Wax.Spike.4.22.21~Entry+(1|Rep),REML=TRUE,K21)
anova(ANOVA_spike_K21)
summary(ANOVA_spike_K21)


#lsmeans
K21spikelsmeansbygeno<-lsmeans(ANOVA_spike_K21,"Entry")
#Warning: Cannot use mode = "kenward-roger" because *pbkrtest* package is not installed
print(K21spikelsmeansbygeno)

#write output
setwd("D:/wax/Field2021/GHpop-lsmeans")
write.csv(K21spikelsmeansbygeno,file="K21-GH-wax-spike-lsmeans.csv")




#### K21 WAX LEAF BOTTOM

#ANOVA
ANOVA_leafB_K21<-lmer(Wax.Leaf.Bottom.4.22.21~Entry+(1|Rep),REML=TRUE,K21)
anova(ANOVA_leafB_K21)
summary(ANOVA_leafB_K21)

#lsmeans
#(LAleafB.rg <- ref.grid(ANOVA_leafB))
K21leafBmeansbygeno<-lsmeans(ANOVA_leafB_K21,"Entry")
print(K21leafBmeansbygeno)

#write output
write.csv(K21leafBmeansbygeno,file="K21-HG-wax-leaf-bottom-lsmeans.csv")
