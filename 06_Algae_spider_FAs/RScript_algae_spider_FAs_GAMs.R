###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Fig.S7 in "Fatty acid biomarkers reveal landscape influences on               ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 10 January 2021 
# Updated: 7 June 2025 

## This R script performs analyses and creates figures relating to the relationship between stream algal accrual and
## riparian spider fatty acid composition. It performs Generalized Additive Model (GAM) analyses to assess the potential
## for non-linear "subsidy-stress" relationships between algal production and realized ecological connectivity

## Metadata for "Burdon_et_al_Fig.S7_multivariate_data.csv" (Numbers = column)
## 1.     Genus_family: Family or taxonomic groups used for arachnids	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Taxa_code: Short name for arachnid family or taxonomic group	
## 4.     Type_1:	Site type (forest, unbuffered, buffered)
## 5.     Type_2:	Site type (paired or reference)
## 6.     Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 7.     Site_no: Site number	
## 8.     Site_name: Site name	
## 9.     Group: Arachnida	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest_level: Lowest taxonomic level used	
## 12.    Tribe:	Used for ground beetles not included in this dataset
## 13.    Ratio_DHA_LNA: Ratio of DHA to LNA
## 14.    Ratio_EPA_LNA: Ratio of EPA to LNA
## 15.    Ratio_EPA_ALA: Ratio of EPA to ALA
## 16.    EPA_perc: %EPA, Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 17.    DHA_perc: %DHA, Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 18.    LNA_perc: %LNA, Linoleic acid (LNA), an omega-6 fatty acid
## 19.    ALA_perc: #ALA, Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 20.    chl_a_m2_day:	Algal accrual on tiles
## 21.    Riparian_PC1:	Riparian condition PC1 (36.7%)
## 22     Catchment_PC1: Catchment impact PC1 (42.3%)	

### Load R packages for data management
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(readr)

### Load R packages for community analysis
library(vegan)

### For data visualization
library(gplots)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(Cairo)
#library(rJava)
#library(venneuler)
library(RColorBrewer)
#install.packages('devtools')
#library(devtools)
#install_github('fawda123/ggord')
library(ggord)

## Testing effects
library(lme4)
library(lmerTest)
library(pbkrtest)
library(lsmeans)
library(MCMCglmm)
library(nlme)
library(phia)
library(lmtest)
library(blme)
library(languageR)

## Loads mgcv package for fitting GAMs
library(mgcv)

### Pairwise test differences in FA composition using pairwise PERMANOVA
#install.packages('devtools')
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

##Loads package used for Logit transformation (% data)
library(car)

##Loads vegan package for transformation
library(vegan)

##Loads packages for correlation check
library(gclus)
library(pheatmap)

## Check if needed
#install.packages("packfor", repos="http://R-Forge.R-project.org", type="source")
#library(packfor)

#############################################################################################
##                                                                                         ##
##                                       LOAD DATA                                         ##        
##                                                                                         ##    
#############################################################################################

## Load data for analyses
Spider_FA_data_new <- read_csv("Burdon_et_al_Fig.S7_multivariate_data.csv")

########################################################################################
##                                                                                    ##
##                      UNIVARIATE RESPONSES FOR KEY VARIABLES                        ##
##                                                                                    ## 
########################################################################################
########################################################################################
##                                                                                    ##
##                           SPIDER: EPA vs SITE vs ALGAE                             ##
##                                                                                    ## 
########################################################################################

## Create new data frame for transformations
Spider_FA_data_log <- Spider_FA_data_new

## Transform percentages using asin sqrt transformation
#Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- Spider_FA_data_log[,-12]

## Remove NAs (256 -> 240 samples)
Spider_FA_data_algae <- na.omit(Spider_FA_data_log)

## Rename data
Spider_FA_data_log[which(Spider_FA_data_new$Type_1=="FOR"),4] <- "3_FOR"
Spider_FA_data_log[which(Spider_FA_data_new$Type_1=="FBF"),4] <- "2_FBF"
Spider_FA_data_log[which(Spider_FA_data_new$Type_1=="UBF"),4] <- "1_UBF"


## Create transformed response variable EPA
Spider_FA_data_log$EPA_perc_asqrt <- asin(sqrt(Spider_FA_data_log$EPA_perc/100))

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(asin(sqrt(Spider_FA_data_log$EPA_perc/100))~Spider_FA_data_log$Type_1,
        col = c("yellowgreen","olivedrab","forestgreen"),
        names = c("No buffer","Forested buffer","Forest"),
        ylab = "% EPA",
        xlab = "Sites type")
boxplot(asin(sqrt(Spider_FA_data_log$EPA_perc/100))~Spider_FA_data_log$Country)

head(Spider_FA_data_log)

## ANOVA: 3 site types
Model_1 <- aov(EPA_perc_asqrt ~ Country/Site_block + Taxa_code*as.factor(Type_1), data = Spider_FA_data_log)
Model_2 <- aov(EPA_perc_asqrt ~ Country/Site_block + Taxa_code + as.factor(Type_1), data = Spider_FA_data_log)
summary(Model_1)
summary(Model_2)
BIC(Model_1)
BIC(Model_2)

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_2b <- aov(EPA_perc_asqrt ~ Country + Taxa_code + as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_2b, ~ Type_1))
test(L.S, adjust = "tukey")

## Need to check further
L.S <- pairs(lsmeans(Model_2b, ~ Country))
test(L.S, adjust = "tukey")

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$EPA_perc_asqrt <- asin(sqrt(Spider_FA_data_algae$EPA_perc/100))
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

## Note: N=240/30

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m2 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 3, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m3 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 4, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m4 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 5, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m5 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 6, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m6 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 7, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m7 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m8 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 9, bs = "cs"), method="REML", data=Spider_FA_data_algae)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m5$gcv.ubre-m6$gcv.ubre ## 0.1649691  M6 k=7 is the best model for log_Algae

summary(m6) 
gam.check(m6, k.sample=5000, k.rep=200, pch=19, cex=.5) ## Looks ok

## Note that "cs" is the cubic spline - bs = "cs" specifies a shrinkage version of "cr".
## These have a cubic spline basis defined by a modest sized set of knots spread evenly through the covariate values.
## They are penalized by the conventional intergrated square second derivative cubic spline penalty.
## More info from: smooth.construct.cr.smooth.spec  - s(x,bs="cs") specifies a penalized cubic regression spline 
## which has had its penalty modified to shrink towards zero at high enough smoothing parameters 
## (as the smoothing parameter goes to infinity a normal cubic spline tends to a straight line.)
## For details see cubic.regression.spline and e.g. Wood (2006a). See also:
## Wood S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition). Chapman and Hall/CRC Press.
## Useful info: https://stats.stackexchange.com/questions/190172/how-i-can-interpret-gam-results
## And more useful info: https://stats.stackexchange.com/questions/359568/choosing-k-in-mgcvs-gam

## Get GAM model for %EPA vs Algae
GAM_EPA_perc_m6 <- gam(EPA_perc_asqrt ~ s(log_Algae, k = 7, bs = "cs"), method="REML", data=Spider_FA_data_algae)

## View relationship
plot(GAM_EPA_perc_m6, ylab = "Arachnid %EPA",
     xlab = expression("log Algae (Chl-a g m"^-2*" day"^-1*")"))
grid.text(expression(paste("s(log Algae,3.38)",sep="")), 
          x = unit(0.15, "npc"), y = unit(0.83, "npc"),just="left", gp=gpar(fontsize=14, col="black"))
grid.text(expression(paste(italic(F)," = 14.9, ",italic(P)," < ",0.001,sep="")), 
          x = unit(0.15, "npc"), just="left", y = unit(0.79, "npc"), gp=gpar(fontsize=14, col="black"))
grid.text(expression(paste("Deviance expl. = 28.4%",sep="")), 
          x = unit(0.15, "npc"), y = unit(0.75, "npc"),just="left", gp=gpar(fontsize=14))

## Linear model
lm1 <- lm(EPA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## RESULT is significant

BIC(lm1)-BIC(m6) ## GAM is an improvement on linear model (note that non-independencies are not accounted for)
AIC(lm1)-AIC(m6)
AICcmodavg::AICc(lm1)-AICcmodavg::AICc(m6)

## Function to calculation VIF in lmer
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam 
  v }

lm3 <- lmer(EPA_perc_asqrt ~ log_Algae*Country + Type_1 + (1|Site_block) + (1|Genus_family), REML=T, data=Spider_FA_data_algae)

summary(lm3)
sjPlot::tab_model(lm3, digits = 5)
car::Anova(lm3,"III")
BIC(lm3)
vif.mer(lm3) ## See variance inflation

lm4 <- lmer(EPA_perc_asqrt ~ log_Algae + Country + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm4)
sjPlot::tab_model(lm4, digits = 5)
car::Anova(lm4,"II")
BIC(lm4)
anova(lm3,lm4) ## No difference with interaction
vif.mer(lm4) ## See variance inflation (Romania and Sweden > 2)

## See: Zuur, A. F., et al. (2010). A protocol for data exploration to avoid common statistical problems.
## Methods in Ecology and Evolution 1(1): 3-14.
## "High, or even moderate, collinearity is especially
## problematic when ecological signals are weak. In that case,
## even a VIF of 2 may cause nonsignificant parameter estimates,
## compared to the situation without collinearity."

## Drop Country > VIF 2
lm2 <- lmer(EPA_perc_asqrt ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
anova(lm4,lm2) ## Country makes a difference, but see VIF > 2
vif.mer(lm2) ## See variance inflation < 2

BIC(m6)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m6)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m6)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(EPA_perc_asqrt ~ log_Algae, 
     ylab = "%EPA",
     xlab = "log Algae",
     main = "SPIDER %EPA",
     data=Spider_FA_data_algae)
abline(lm(EPA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae), col="red")
summary(lm(EPA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae))

##############################################################################
##              LOOK AT STABILISING INFLUENCE OF BUFFER                     ##
##############################################################################

plot(EPA_perc_asqrt ~ Riparian_PC1, 
     xlab = "Riparian_PC1",
     ylab = "%EPA",
     main = "SPIDER %EPA",
     data=Spider_FA_data_algae)
abline(lm(EPA_perc_asqrt ~ Riparian_PC1, data=Spider_FA_data_algae), col="red")
summary(lm(EPA_perc_asqrt ~ Riparian_PC1, data=Spider_FA_data_algae))

Linear_M1 <- lm(EPA_perc_asqrt ~ Riparian_PC1, data=Spider_FA_data_algae)
Linear_M1_residual <- resid(Linear_M1)
Residual_deviation <- Linear_M1_residual*Linear_M1_residual 

plot(Residual_deviation ~ Riparian_PC1, 
     xlab = "Riparian_PC1",
     ylab = "Residual deviation",
     main = "SPIDER %EPA",
     data=Spider_FA_data_algae)
abline(lm(Residual_deviation  ~ Riparian_PC1, data=Spider_FA_data_algae), col="red")
summary(lm(Residual_deviation  ~ Riparian_PC1, data=Spider_FA_data_algae)) ## Evidence that Riparian PC1 stabilizes EPA%

head(Spider_FA_data_algae)
Model <- lmer(Residual_deviation ~ Riparian_PC1 + Type_1 + (1|Site_block) + (1|Genus_family), data=Spider_FA_data_algae) ## Evidence that Riparian PC1 stabilizes EPA%

summary(Model)
sjPlot::tab_model(Model , digits = 5)
car::Anova(Model ,"II")
BIC(Model)
vif.mer(Model) ## See variance inflation < 2

Modelb <- lmer(Residual_deviation ~ Riparian_PC1*Country + Type_1 + (1|Site_block) + (1|Genus_family), data=Spider_FA_data_algae) ## Evidence that Riparian PC1 stabilizes EPA%

summary(Modelb)
sjPlot::tab_model(Modelb, digits = 5)
car::Anova(Modelb,"III")
BIC(Modelb)
vif.mer(Modelb) ## See variance inflation < 2 (Drop Country)

Modelc <- lmer(Residual_deviation ~ Riparian_PC1 + Country + Type_1 + (1|Site_block) + (1|Genus_family), data=Spider_FA_data_algae) ## Evidence that Riparian PC1 stabilizes EPA%

summary(Modelc)
sjPlot::tab_model(Modelc, digits = 5)
car::Anova(Modelc,"II")
BIC(Modelc)
vif.mer(Modelc) ## See variance inflation < 2 (Drop Country)

anova(Modelb,Modelc) ## See interaction

########################################################################################
##                                                                                    ##
##                           SPIDER: DHA vs SITE vs ALGAE                             ##
##                                                                                    ## 
########################################################################################

row.names(Spider_FA_data_log) <- NULL

## Create transformed response variable EPA
Spider_FA_data_log$DHA_perc_asqrt <- asin(sqrt(Spider_FA_data_log$DHA_perc/100))

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(asin(sqrt(Spider_FA_data_log$DHA_perc/100))~Spider_FA_data_log$Type_1)
boxplot(asin(sqrt(Spider_FA_data_log$DHA_perc/100))~Spider_FA_data_log$Country)

## ANOVA - three extreme values forcing difference between site types
Model_1 <- aov(DHA_perc_asqrt ~ as.factor(Type_1) + Taxa_code + Country/Site_block, data = Spider_FA_data_log[-c(132,187,239),])
summary(Model_1)
BIC(Model_1)
plot(Model_1)

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_1b <- glm(DHA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## Write post-hoc test output GLM
Model_1b <- glm(DHA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$DHA_perc_asqrt <- asin(sqrt(Spider_FA_data_algae$DHA_perc/100))
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m2 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m3 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m4 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m5 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m6 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m7 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m8 <- gam(DHA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m2$gcv.ubre-m1$gcv.ubre ## M1 k=8 no "one" best model for log_Algae

summary(m1) 
gam.check(m1)
plot(m1, ylab = "s(log Algae,1)",
     xlab = "log Algae",
     main = " SPIDER %DHA")

lm1 <- lm(DHA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## Note that the relationship is non-significant
BIC(lm1)-BIC(m1) ## GAM is NOT an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(DHA_perc_asqrt ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m1)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m1)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m1)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(DHA_perc_asqrt ~ log_Algae, 
     ylab = "%DHA",
     xlab = "log Algae",
     main = "SPIDER %DHA",
     data=Spider_FA_data_algae)
abline(lm(DHA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae), col="red")

########################################################################################
##                                                                                    ##
##                         SPIDER: ALA vs SITE vs ALGAE                               ##
##                                                                                    ## 
########################################################################################

## Create transformed response variable EPA
Spider_FA_data_log$ALA_perc_asqrt <- asin(sqrt(Spider_FA_data_log$ALA_perc/100))

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(asin(sqrt(Spider_FA_data_log$ALA_perc/100))~Spider_FA_data_log$Type_1)
boxplot(asin(sqrt(Spider_FA_data_log$ALA_perc/100))~Spider_FA_data_log$Country)

## ANOVA
Model_1 <- aov(ALA_perc_asqrt ~ as.factor(Type_1) + Taxa_code + Country/Site_block, data = Spider_FA_data_log)
summary(Model_1)
BIC(Model_1)
#plot(Model_1)

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_1b <- glm(ALA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## Write post-hoc test output GLM
Model_1b <- glm(ALA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$ALA_perc_asqrt <- asin(sqrt(Spider_FA_data_algae$ALA_perc/100))
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m2 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=3, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m3 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=4, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m4 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=5, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m5 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=6, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m6 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=7, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m7 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m8 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=9, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m9 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=10, bs = "cs"), method="REML",data=Spider_FA_data_algae)
m10 <- gam(ALA_perc_asqrt ~ s(log_Algae, k=11, bs = "cs"), method="REML",data=Spider_FA_data_algae)


m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  
m9$gcv.ubre[1]  
m10$gcv.ubre[1]  

m7$gcv.ubre-m5$gcv.ubre ## dREML 0.226359  M5 k=6 is the best model for log_Algae

summary(m9) 
gam.check(m9)
plot(m4, ylab = "s(log Algae,1)",
     xlab = "log Algae",
     main = " SPIDER %ALA")

lm1 <- lm(ALA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## Note that the relationship is non-significant

BIC(lm1)-BIC(m4) ## GAM is an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(ALA_perc_asqrt ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m5)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m5)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m5)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(ALA_perc_asqrt ~ log_Algae, 
     ylab = "%ALA",
     xlab = "log Algae",
     main = "SPIDER %ALA",
     data=Spider_FA_data_algae)
abline(lm(ALA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae), col="red")

########################################################################################
##                                                                                    ##
##                           SPIDER: LNA vs SITE vs ALGAE                             ##
##                                                                                    ## 
########################################################################################

## Create transformed response variable LNA
Spider_FA_data_log$LNA_perc_asqrt <- asin(sqrt(Spider_FA_data_log$LNA_perc/100))

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(asin(sqrt(Spider_FA_data_log$LNA_perc/100))~Spider_FA_data_log$Type_1)
boxplot(asin(sqrt(Spider_FA_data_log$LNA_perc/100))~Spider_FA_data_log$Country)

## ANOVA
Model_1 <- aov(LNA_perc_asqrt ~ as.factor(Type_1) + Taxa_code + Country/Site_block, data = Spider_FA_data_log)
summary(Model_1)
BIC(Model_1)
#plot(Model_1)

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_1b <- glm(LNA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## Write post-hoc test output GLM
Model_1b <- glm(LNA_perc_asqrt ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$LNA_perc_asqrt <- asin(sqrt(Spider_FA_data_algae$LNA_perc/100))
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m2 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=3, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m3 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=4, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m4 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=5, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m5 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=6, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m6 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=7, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m7 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m8 <- gam(LNA_perc_asqrt ~ s(log_Algae, k=9, bs = "cs"), method="REML", data=Spider_FA_data_algae)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m7$gcv.ubre-m6$gcv.ubre ## dgcv.ubre 0.1309197 M6 k=7 is the best model for log_Algae

summary(m6) 
gam.check(m6)
plot(m5, ylab = "s(log Algae,2.42)",
     xlab = "log Algae",
     main = " SPIDER %LIN")

lm1 <- lm(LNA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## Note that the relationship is non-significant

BIC(lm1)-BIC(m6) ## GAM is an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(LNA_perc_asqrt ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m6)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m6)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m6)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(LNA_perc_asqrt ~ log_Algae, 
     ylab = "%LIN",
     xlab = "log Algae",
     main = "SPIDER %LIN",
     data=Spider_FA_data_algae)
abline(lm(LNA_perc_asqrt ~ log_Algae, data=Spider_FA_data_algae), col="red")

########################################################################################
##                                                                                    ##
##                    SPIDER: RATIO OF DHA to LIN vs SITE vs ALGAE                    ##
##                                                                                    ## 
########################################################################################

## CHECK RATIO OF DHA to LNA
Spider_FA_data_log$Ratio_DHA_LNA

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(log(Spider_FA_data_log$Ratio_DHA_LNA)~Spider_FA_data_log$Type_1)
boxplot(log(Spider_FA_data_log$Ratio_DHA_LNA)~Spider_FA_data_log$Country)

## ANOVA
Model_1 <- aov(log(Ratio_DHA_LNA) ~ Type_1 + Taxa_code + Country/Site_block, data = Spider_FA_data_log)
summary(Model_1)
BIC(Model_1)
#plot(Model_1)

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_1b <- aov(log(Ratio_DHA_LNA) ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## Write post-hoc test output GLM
Model_1b <- glm(log(Ratio_DHA_LNA) ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m2 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=3, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m3 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=4, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m4 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=5, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m5 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=6, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m6 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=7, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m7 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m8 <- gam(log(Ratio_DHA_LNA) ~ s(log_Algae, k=9, bs = "cs"), method="REML", data=Spider_FA_data_algae)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m8$gcv.ubre-m1$gcv.ubre ## M1 k=8 is the best model for log_Algae

summary(m1) 
gam.check(m1)
plot(m1, ylab = "s(log Algae,0.80)",
     xlab = "log Algae",
     main = " SPIDER DHA:LNA")

lm1 <- lm(log(Ratio_DHA_LNA) ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## Note that the + relationship is non-significant

BIC(lm1)-BIC(m4) ## GAM is NOT an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(log(Ratio_DHA_LNA) ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m1)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m1)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m1)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(log(Ratio_DHA_LNA) ~ log_Algae, 
     ylab = "DHA:LIN",
     xlab = "log Algae",
     main = "SPIDER DHA:LIN",
     data=Spider_FA_data_algae)
abline(lm(log(Ratio_DHA_LNA) ~ log_Algae, data=Spider_FA_data_algae), col="red")

########################################################################################
##                                                                                    ##
##                   SPIDER: RATIO OF EPA to ALA vs SITE vs ALGAE                     ##
##                                                                                    ## 
########################################################################################

## Here the hypothesis is that ALA is a potential precursor of EPA, thus there should be
## an inverse (opposite) of the ratio for EPA:ALA as connectivity increases

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(Spider_FA_data_log$Ratio_EPA_ALA~Spider_FA_data_log$Type_1,
        col = c("yellowgreen","olivedrab","forestgreen"),
        names = c("No buffer","Forested buffer","Forest"),
        ylab = "EPA:ALA",
        xlab = "Site Type")
boxplot(Spider_FA_data_log$Ratio_EPA_ALA~Spider_FA_data_log$Country)

head(Spider_FA_data_log)

## ANOVA
Model_1 <- aov(Ratio_EPA_ALA ~ Country/Site_block + Mode + Taxa_code + Type_1, data = Spider_FA_data_log)
summary(Model_1)
BIC(Model_1)

## Calculate mean and median values to interpret results of above
Spider_FA_data_log %>% 
  group_by(Type_1) %>% summarize(mean = mean(Ratio_EPA_ALA, na.rm = TRUE),
                                 sd = sd(Ratio_EPA_ALA, na.rm = TRUE),
                                 median = median(Ratio_EPA_ALA, na.rm = TRUE))

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
Model_1b <- aov(Ratio_EPA_ALA ~ as.factor(Type_1), data = Spider_FA_data_log)
L.S <- pairs(lsmeans(Model_1b, ~ Type_1))
test(L.S, adjust = "tukey")

## ANOVA - site pairs to better understand the difference indicated above
Model_1 <- aov(Ratio_EPA_ALA ~ Country/Site_block + Mode + Taxa_code + Type_1, data = Spider_FA_data_log[which(Spider_FA_data_log$Type_2=="PAIR"),])
summary(Model_1)
BIC(Model_1)

## ALGAE

## Create Algae predictor variable
Spider_FA_data_algae$log_Algae <- log1p(Spider_FA_data_algae$chl_a_m2_day)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m2 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=3, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m3 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=4, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m4 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=5, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m5 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=6, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m6 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=7, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m7 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae)
m8 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=9, bs = "cs"), method="REML", data=Spider_FA_data_algae)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m6$gcv.ubre-m2$gcv.ubre ## M2 k=3 is the best model for log_Algae

summary(m2) 
gam.check(m2)
plot(m2, ylab = "s(log Algae,0.79)",
     xlab = "log Algae",
     main = " SPIDER EPA:ALA")

lm1 <- lm(Ratio_EPA_ALA ~ log_Algae, data=Spider_FA_data_algae)
summary(lm1) ## Note that the + relationship is significant
#plot(lm1)
#resid(lm1)

BIC(lm1)-BIC(m2) ## GAM is an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(log(Ratio_EPA_ALA) ~ log_Algae + Type_1 + (1|Site_block) + (1|Genus_family), REML=F, data=Spider_FA_data_algae)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

## Create object for selected model
GAM_EPA_ALA_m2 <- gam(Ratio_EPA_ALA ~ s(log_Algae, k=3, bs = "cs"), method="REML", data=Spider_FA_data_algae)

BIC(m2)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m2)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m2)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(Ratio_EPA_ALA ~ log_Algae, 
     ylab = "EPA:ALA",
     xlab = "log Algae",
     main = "SPIDER EPA:ALA",
     data=Spider_FA_data_algae)
abline(lm(Ratio_EPA_ALA ~ log_Algae, data=Spider_FA_data_algae), col="red")


########################################################################################
##                                                                                    ##
##                      SPIDER: ALGAE VS TROPHIC CONNECTIVITY                         ##
##                                                                                    ## 
########################################################################################

## Create plot for %EPA and EPA:ALA
png(filename="Burdon_et_al_Fig.S7.png", 
    type="cairo",
    units="in", 
    width=12, 
    height=8, 
    pointsize=20, 
    res=600)

par(mfcol=c(1,2))

plot(GAM_EPA_perc_m6, ylab = "Arachnid %EPA",
     xlab = expression("log Algae (Chl-a g m"^-2*" day"^-1*")"))
grid.text(expression(paste("s(log Algae,3.38)",sep="")), 
          x = unit(0.13, "npc"), y = unit(0.80, "npc"),just="left", gp=gpar(fontsize=12, col="black"))
grid.text(expression(paste(italic(F)," = 14.9, ",italic(P)," < ",0.001,sep="")), 
          x = unit(0.13, "npc"), just="left", y = unit(0.76, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text(expression(paste("Deviance expl. = 28.4%",sep="")), 
          x = unit(0.13, "npc"), y = unit(0.725, "npc"),just="left", gp=gpar(fontsize=12))

plot(GAM_EPA_ALA_m2, ylab = "Arachnid EPA:ALA",
     xlab = expression("log Algae (Chl-a g m"^-2*" day"^-1*")"))
grid.text(expression(paste("s(log Algae,1.79)",sep="")), 
          x = unit(0.63, "npc"), y = unit(0.80, "npc"),just="left", gp=gpar(fontsize=12, col="black"))
grid.text(expression(paste(italic(F)," = 20.5, ",italic(P)," < ",0.001,sep="")), 
          x = unit(0.63, "npc"), just="left", y = unit(0.76, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text(expression(paste("Deviance expl. = 15.4%",sep="")), 
          x = unit(0.63, "npc"), y = unit(0.725, "npc"),just="left", gp=gpar(fontsize=12))

grid.text(expression(paste("(a)",sep="")), 
          x = unit(0.03, "npc"), y = unit(0.82, "npc"),just="centre", gp=gpar(fontsize=16, fontface="bold"))
grid.text(expression(paste("(b)",sep="")), 
          x = unit(0.53, "npc"), y = unit(0.82, "npc"),just="centre", gp=gpar(fontsize=16, fontface="bold"))

dev.off()

########################################################################################
##                                                                                    ##
##                      SPIDER: ALGAE VS RIPARIAN CONDITION                           ##
##                                                                                    ## 
########################################################################################

## Use site means to avoid confounding spider presence with site algae
Spider_FA_data_algae_mu <- Spider_FA_data_algae %>% 
  group_by(Country,Type_1,Type_2,Site_block,Site_no,Site_name) %>% 
  summarise(Algae=mean(chl_a_m2_day),
            Ratio_EPA_ALA=mean(Ratio_EPA_ALA),
            Riparian_PC1=mean(Riparian_PC1))

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=3, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m2 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=3, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m3 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=4, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m4 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=5, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m5 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=6, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m6 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=7, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m7 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=8, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)
m8 <- gam(log1p(Algae) ~ s(Riparian_PC1, k=9, bs = "cs"), method="REML",  data=Spider_FA_data_algae_mu)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m2$gcv.ubre-m1$gcv.ubre ## M1 k=3 is the best model for log_Algae

summary(m1) 
gam.check(m1)
plot(m1, ylab = "s(Riparian_PC1,0.000035)",
     xlab = "Riparian_PC1",
     main = "SPIDER log Algae")

lm1 <- lm(log1p(Algae) ~ Riparian_PC1, data=Spider_FA_data_algae_mu)
summary(lm1) ## Note that the + relationship is non-significant
BIC(lm1)-BIC(m1) ## GAM is NOT an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(log1p(Algae) ~ Riparian_PC1 + Type_1 + (1|Site_block), REML=F, data=Spider_FA_data_algae_mu)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m1)-BIC(lm2) ## GAM is not an improvement on mixed model
AIC(m1)-AIC(lm2) ## GAM is not an improvement on mixed model
AICcmodavg::AICc(m1)-AICcmodavg::AICc(lm2) ## GAM is not an improvement on mixed model

plot(log1p(Algae) ~ Riparian_PC1, 
     xlab = "Riparian_PC1",
     ylab = "log Algae",
     main = "SPIDER log Algae",
     data=Spider_FA_data_algae_mu)
abline(lm(log1p(Algae) ~ Riparian_PC1, data=Spider_FA_data_algae_mu), col="red", lty=2)

########################################################################################
##                                                                                    ##
##                 SPIDER: RATIO OF EPA to ALA VS RIPARIAN CONDITION                  ##
##                                                                                    ## 
########################################################################################

## SITE TYPE (FOREST vs FBF vs UBF)
boxplot(Spider_FA_data_log$Riparian_PC1~Spider_FA_data_log$Type_1)

### Assess potential for non-linear relationship between Algae production and trophic connectivity
m1 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m2 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=3, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m3 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=4, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m4 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=5, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m5 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=6, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m6 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=7, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m7 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=8, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)
m8 <- gam(Ratio_EPA_ALA ~ s(Riparian_PC1, k=9, bs = "cs"), method="REML", data=Spider_FA_data_algae_mu)

m1$gcv.ubre[1] 
m2$gcv.ubre[1] 
m3$gcv.ubre[1] 
m4$gcv.ubre[1]  
m5$gcv.ubre[1]  
m6$gcv.ubre[1] 
m7$gcv.ubre[1]  
m8$gcv.ubre[1]  

m6$gcv.ubre-m1$gcv.ubre ## delta = 0.2713504 M1 k=8 is the best model for log_Algae

summary(m1) 
gam.check(m1)
plot(m1, ylab = "s(Riparian_PC1,3.78)",
     xlab = "Riparian_PC1",
     main = "SPIDER EPA:ALA")

lm1 <- lm(Ratio_EPA_ALA ~ Riparian_PC1, data=Spider_FA_data_algae_mu)
summary(lm1) ## Note that the -ve relationship is non-significant

BIC(lm1)-BIC(m1) ## GAM is NOT an improvement on linear model (note that non-independencies are not accounted for)

## Drop Country > VIF 2
lm2 <- lmer(Ratio_EPA_ALA ~ Riparian_PC1 + Type_1 + (1|Site_block), REML=F, data=Spider_FA_data_algae_mu)

summary(lm2)
sjPlot::tab_model(lm2, digits = 5)
car::Anova(lm2,"II")
BIC(lm2)
vif.mer(lm2) ## See variance inflation < 2

BIC(m1)-BIC(lm2)## Mixed effects improvement on GAM  

plot(Ratio_EPA_ALA ~ Riparian_PC1, 
     xlab = "Riparian_PC1",
     ylab = "EPA:ALA",
     main = "SPIDER EPA:ALA",
     data=Spider_FA_data_algae_mu)
abline(lm(Ratio_EPA_ALA ~ Riparian_PC1, data=Spider_FA_data_algae_mu), col="red", lty=2)