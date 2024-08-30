## Entry means calculation for wax ratings
# HilliardxGA-13LE6 

#install.packages("lmerTest")

library("lme4")
library("lmerTest")
library("emmeans")

#### import phenotype combine single dataset ####

### all wax trts R20 LA pop
R20waxLA<-read.table("C:/Users/Daniela/Documents/Wax/wax_mini/Field2020/R20_LA_allwax_wparents.csv",sep=",",head=TRUE,
                     colClasses=c("numeric","character","numeric","numeric","numeric","numeric"), na.strings=".")
head(R20waxLA, n=5)
tail(R20waxLA, n=5)
colnames(R20waxLA)
#rename cols
R20waxLA$Wax.Spike <- R20waxLA$WaxHeads
R20waxLA$Wax.Leaf <- R20waxLA$WaxLeafBOTTOM
R20waxLA <- subset (R20waxLA, select = c(-WaxHeads, -WaxLeafBOTTOM, -WaxLeafTOP, -WaxLeafComp))
R20waxLA$Env <- "RAL20"
colnames(R20waxLA)


### all wax trts R21 LA pop
R21waxLA<-read.table("C:/Users/Daniela/Documents/Wax/wax_mini/Field2021/R21_T26-T34_LA_alltrts_wparents.csv",sep=",",head=TRUE,
                     colClasses=c("numeric","character","numeric","numeric","character"), na.strings=".")
head(R21waxLA, n=5)
tail(R21waxLA, n=5)
colnames(R21waxLA)
R21waxLA$Wax.Spike <- R21waxLA$Wax.Spike.4.26.21
R21waxLA$Wax.Leaf <- R21waxLA$Wax.Leaf.Bottom.4.22.21
R21waxLA <- subset (R21waxLA, select = c(-Wax.Spike.4.26.21, -Wax.Leaf.Bottom.4.22.21, -flowering))
R21waxLA$Env <- "RAL21"
colnames(R21waxLA)

## combine into single dataset
LAcom <- rbind(R20waxLA, R21waxLA)
colnames(LAcom)
head(LAcom)
tail(LAcom)

#make sure rep is a factor variable
LAcom$Rep <- as.factor(LAcom$Rep)
LAcom$Env <- as.factor(LAcom$Env)
class(LAcom$Rep)
class(LAcom$Env)


#####  LA spike combined env model #####

#using same model as HG pop
require(lme4)
mm_spike <- lmer(Wax.Spike ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, LAcom)
summary(mm_spike) 
isSingular(mm_spike) #SINGULAR
anova(mm_spike)
#overfit model -- use simpler model without Env effects

#model without Env effects
require(lme4)
mm_spike1 <- lmer(Wax.Spike ~ 1 + Entry + (1|Rep), REML=TRUE, LAcom)
summary(mm_spike1) 
isSingular(mm_spike1) #NOT SINGULAR
anova(mm_spike1)
# -- use this model 

#calculate Entry means
emspikeLA <- emmeans(mm_spike1,"Entry")
print(emspikeLA)

#write output
setwd("~/Wax/LAfinal/means/")
write.csv(emspikeLA,file="spike-emmeans-across-envs.csv")


#####  LA leaf combined env model #####

#using same model as HG pop
require(lme4)
mm_leaf <- lmer(Wax.Leaf ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, LAcom)
summary(mm_leaf) 
isSingular(mm_leaf) #SINGULAR
anova(mm_leaf)
#overfit model, remove Env effects

#model without Env effects
require(lme4)
mm_leaf1 <- lmer(Wax.Leaf ~ 1 + Entry + (1|Rep), REML=TRUE, LAcom)
summary(mm_leaf1) 
isSingular(mm_leaf1) #NOT SINGULAR
anova(mm_leaf1)
# -- use this model 

#calculate Entry means
emleafLA <- emmeans(mm_leaf1,"Entry")
print(emleafLA)

#write output
setwd("~/Wax/LAfinal/means/")
write.csv(emleafLA,file="leaf-emmeans-across-envs.csv")


#### R21 WAX LEAF #### 

#ANOVA
ANOVA_leafB<-lmer(Wax.Leaf~Entry+(1|Rep),REML=TRUE,R21waxLA)
anova(ANOVA_leafB)
summary(ANOVA_leafB)

#lsmeans
leafBmeansbygeno<-lsmeans(ANOVA_leafB,"Entry")
print(leafBmeansbygeno)

#write output
write.csv(leafBmeansbygeno,file="R21-LA-wax-leaf-bottom-lsmeans.csv")


#### R21 WAX SPIKE #### 

#ANOVA
ANOVA_spike<-lmer(Wax.Spike~Entry+(1|Rep),REML=TRUE,R21waxLA)
anova(ANOVA_spike)
summary(ANOVA_spike)

#lsmeans
spikemeansbygeno<-lsmeans(ANOVA_spike,"Entry")
print(spikemeansbygeno)

#write output
write.csv(spikemeansbygeno,file="R21-LA-wax-spike-lsmeans.csv")


#### R20 WAX LEAF #### 

#ANOVA
ANOVA_leafBK<-lmer(Wax.Leaf~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_leafBK)
summary(ANOVA_leafBK)

#lsmeans
leafBmeansbygenoK<-lsmeans(ANOVA_leafBK,"Entry")
print(leafBmeansbygenoK)

#write output
write.csv(leafBmeansbygenoK,file="R20-LA-wax-leaf-bottom-lsmeans.csv")


#### R20 WAX SPIKE #### 

#ANOVA
ANOVA_spikeK<-lmer(Wax.Spike~Entry+(1|Rep),REML=TRUE,R20waxLA)
anova(ANOVA_spikeK)
summary(ANOVA_spikeK)

#lsmeans
spikemeansbygenoK<-lsmeans(ANOVA_spikeK,"Entry")
print(spikemeansbygenoK)

#write output
write.csv(spikemeansbygenoK,file="R20-LA-wax-spike-lsmeans.csv")

