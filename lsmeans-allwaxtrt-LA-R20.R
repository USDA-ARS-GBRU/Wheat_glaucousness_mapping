## Lsmeans calculation for wax leaf types
#Raleigh 2020 LA population

library("lme4")
library("lmerTest")
library("lsmeans")

### all wax trts R20 LA pop
R20waxLA<-read.table("D:/wax/Field2020/R20_LA_allwax_wparents.csv",sep=",",head=TRUE,
                     colClasses=c("numeric","character","numeric","numeric","numeric","numeric"), na.strings=".")
head(R20waxLA, n=5)
tail(R20waxLA, n=5)
colnames(R20waxLA)



#### wAX SPIKE ####

#ANOVA
###Type3 Test of hypothesis of fixed effects
###Note: Use Year and Line as fixed effect
ANOVA_R20LAspike<-lmer(WaxHeads~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_R20LAspike, type="III")
summary(ANOVA_R20LAspike)

#lsmeans
(LAwaxheadsR.rg <- ref.grid(ANOVA_R20LAspike))
R20LAwaxheadslsmeansbygeno<-lsmeans(ANOVA_R20LAspike,"Entry")
print(R20LAwaxheadslsmeansbygeno)

#write output
setwd("D:/wax/Field2020/lsmeans")
write.csv(R20LAwaxheadslsmeansbygeno,file="R20-LA-wax-spike-lsmeans.csv")




#### WAX LEAF BOTTOM #### 

#ANOVA
ANOVA_R20LAleafB<-lmer(WaxLeafBOTTOM~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_R20LAleafB, type="III")
summary(ANOVA_R20LAleafB)

#lsmeans
(LAleafB.rg <- ref.grid(ANOVA_R20LAleafB))
R20LAleafBmeansbygeno<-lsmeans(ANOVA_R20LAleafB,"Entry")
print(R20LAleafBmeansbygeno)

#write output
write.csv(R20LAleafBmeansbygeno,file="R20-LA-wax-leaf-bottom-lsmeans.csv")





#### WAX LEAF TOP ####

#ANOVA
ANOVA_R20LAleafT<-lmer(WaxLeafTOP~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_R20LAleafT, type="III")
summary(ANOVA_R20LAleafT)

#lsmeans
(LAleafT.rg <- ref.grid(ANOVA_R20LAleafT))
R20LAleafTmeansbygeno<-lsmeans(ANOVA_R20LAleafT,"Entry")
print(R20LAleafTmeansbygeno)

#write output
write.csv(R20LAleafTmeansbygeno,file="R20-LA-wax-leaf-top-lsmeans.csv")





#### WAX LEAF COMPOSITE ####

# wax leaf composite = wax leaf bottom + wax leaf top

#ANOVA
ANOVA_R20leafcompLA<-lmer(WaxLeafComp~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_R20leafcompLA, type="III")
summary(ANOVA_R20leafcompLA)

#lsmeans
(LAwaxleafcompR.rg <- ref.grid(ANOVA_R20leafcompLA))
R20LAwaxleafcomplsmeansbygeno<-lsmeans(ANOVA_R20leafcompLA,"Entry")
print(R20LAwaxleafcomplsmeansbygeno)

#write output
write.csv(R20LAwaxleafcomplsmeansbygeno,file="R20-LA-wax-leaf-comp-lsmeans.csv")




