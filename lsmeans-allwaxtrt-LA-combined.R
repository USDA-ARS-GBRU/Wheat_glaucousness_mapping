## Lsmeans calculation for wax leaf types
#Raleigh 2020 LA population

library("lme4")
library("lsmeans")

### all wax trts R20 LA pop
R21waxLA<-read.table("D:/wax/Field2021/R21_T26-T34_LA_alltrts_wparents.csv",sep=",",head=TRUE,
                     colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
head(R21waxLA, n=5)
tail(R21waxLA, n=5)
colnames(R21waxLA)

#Julian Date

#set HD col to "Date" class
d <- as.Date(R21waxLA$flowering, format="%m/%d/%Y") #raw formatted in month-day-year
d #default format: year-month-day

#convert "Date" class to "numeric" Julian date 
JulianDay <- julian(d, origin = as.Date("2020-01-01")) #set origin to first of the year from data
JulianDay
class(JulianDay)


#add Julian date col to dataframe
R21waxLA$JulianHD <- JulianDay
head(R21waxLA)

colnames(R21waxLA)

#### wAX SPIKE ####

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
ANOVA_spike<-lmer(Wax.Spike.4.26.21~Entry+(1|Rep),REML=TRUE,R21waxLA)
anova(ANOVA_spike, type="III")
summary(ANOVA_spike)

#lsmeans
(LAwaxheadsR.rg <- ref.grid(ANOVA_spike))
LAspikelsmeansbygeno<-lsmeans(ANOVA_spike,"Entry")
print(LAspikelsmeansbygeno)

#write output
setwd("D:/wax/Field2021/lsmeans")
write.csv(LAspikelsmeansbygeno,file="R21-LA-wax-spike-lsmeans.csv")




#### WAX LEAF BOTTOM #### 

#ANOVA
ANOVA_leafB<-lmer(Wax.Leaf.Bottom.4.22.21~Entry+(1|Rep),REML=TRUE,R21waxLA)
anova(ANOVA_leafB, type="III")
summary(ANOVA_leafB)

#lsmeans
(LAleafB.rg <- ref.grid(ANOVA_leafB))
LAleafBmeansbygeno<-lsmeans(ANOVA_leafB,"Entry")
print(LAleafBmeansbygeno)

#write output
write.csv(LAleafBmeansbygeno,file="R21-LA-wax-leaf-bottom-lsmeans.csv")





#### HEADING DATE - JULIAN DAYS ####

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
ANOVA_HD<-lmer(JulianHD~Entry+(1|Rep),REML=TRUE,R21waxLA)
anova(ANOVA_HD, type="III")
summary(ANOVA_HD)

#lsmeans
(LAHDR.rg <- ref.grid(ANOVA_HD))
LAHDlsmeansbygeno<-lsmeans(ANOVA_HD,"Entry")
print(LAHDlsmeansbygeno)

#write output
write.csv(LAHDlsmeansbygeno,file="R21-LA-HDjulian-lsmeans.csv")



