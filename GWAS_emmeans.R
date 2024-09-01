library("lme4")
library("lmerTest")
library("emmeans")

##### IMPORT DATA #####

#Raleigh 2020 environment
R20<-read.table("C:/Users/Daniela/Documents/Wax/GAPIT/Raleigh2020_wax.csv",sep=",",head=TRUE,
                colClasses=c("character","numeric","numeric","numeric"), na.strings=".")
head(R20, n=5)
tail(R20, n=5)
colnames(R20)
#clean up colnames
R20$Wax.Spike <- R20$WaxSpike
R20$Wax.Leaf <- R20$WaxLeafBOTTOM
R20 <- subset (R20, select = c(-WaxSpike, -WaxLeafBOTTOM))
R20$Env <- "RAL20"
colnames(R20)

#Raleigh 2021 environment
R21<-read.table("C:/Users/Daniela/Documents/Wax/GAPIT/Raleigh2021_wax.csv",sep=",",head=TRUE,
                colClasses=c("character","numeric","numeric","numeric"), na.strings=" ")
head(R21, n=5)
tail(R21, n=5)
colnames(R21)
#clean up colnames
R21$Wax.Spike <- R21$R21Spike
R21$Wax.Leaf <- R21$R21LeafBottom
R21 <- subset (R21, select = c(-R21Spike, -R21LeafBottom))
R21$Env <- "RAL21"
colnames(R21)

##### COMBINE INTO ONE DATASET #####

total <- rbind(R20, R21)
colnames(total)
head(total)
tail(total)

#make sure rep is a factor variable
total$Rep <- as.factor(total$Rep)
#total$Entry <- as.factor(total$Entry)
total$Env <- as.factor(total$Env)

class(total$Entry)
class(total$Rep)
class(total$Wax.Leaf)
class(total$Wax.Spike)
class(total$Env)


##### SPIKE combined MODEL #####

#model including environment effect
require(lme4)
mm_spike7 <- lmer(Wax.Spike ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, total)
summary(mm_spike7) 
isSingular(mm_spike7) #NOT SINGULAR
anova(mm_spike7)

#to estimate environment and rep effects
require(lme4)
mm_spike8 <- lmer(Wax.Spike ~ 1 + Entry + Env + Rep:Env  + (1|Env:Entry), REML=TRUE, total)
isSingular(mm_spike8) #NOT SINGULAR
anova(mm_spike8)
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Entry   449.90  2.5708   175 169.00  8.4931 < 2.2e-16 ***
#Env       4.61  4.6132     1 406.83 15.2403 0.0001108 ***
#Env:Rep   7.37  3.6860     2 344.00 12.1770 7.769e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##### SPIKE combined MEANS #####

emspike <- emmeans(mm_spike7,"Entry")
print(emspike)
head(emspike)
#write output
setwd("C:/Users/Daniela/Documents/Wax/GAPIT/emmeans/")
write.csv(emspike,file="GWAS-spike-emmeans-across-envs.csv")


##### LEAF combined MODEL #####
#leaf glaucousness 

#acrosslocs
require(lme4)
mm_leafg <- lmer(Wax.Leaf ~ 1 + Entry + (1|Env) + (1|Rep:Env)  + (1|Env:Entry), REML=TRUE, total)
isSingular(mm_leafg)
anova(mm_leafg)

#to estimate environment and rep effects
require(lme4)
mm_leaf8 <- lmer(Wax.Leaf ~ 1 + Entry + Env + Rep:Env  + (1|Env:Entry), REML=TRUE, total)
isSingular(mm_leaf8) #NOT SINGULAR
anova(mm_leaf8)
#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Entry   384.53  2.1973   175 169.00  5.3809 < 2.2e-16 ***
#Env      18.00 17.9958     1 348.67 44.0684 1.216e-10 ***
#Env:Rep   3.02  1.5121     2 344.00  3.7029   0.02564 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#simple model no env effect
#require(lme4)
#mm_leaf1 <- lmer(Wax.Leaf ~ 1 + Entry + (1|Rep:Env), REML=TRUE, total)
#isSingular(mm_leaf1)
#anova(mm_leaf1)


##### LEAF combined MEANS #####

emleaf <- emmeans(mm_leafg, "Entry")
print(emleaf)
head(emleaf)
#write output
write.csv(emleaf,file="GWAS-leaf-emmeans-across-envs.csv")

#both traits emmeans combined in a single file "waxtrt_emmeans_GWAS_combenv.csv"