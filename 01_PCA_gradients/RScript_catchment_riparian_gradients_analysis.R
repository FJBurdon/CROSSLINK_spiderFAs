###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Fig.3 in "Fatty acid biomarkers reveal landscape influences on                ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 30 January 2024
# Updated: 7 June 2025

## This R script performs analyses and creates figures relating to the comparison of two environmental gradients
## Environmental gradients constructed using Principal Components Analysis
## The two environmental gradients are based on catchment impacts and riparian condition, respectively
## The script create boxplots by site type and country
## The script also performed linear mixed model analyses to test difference between site types, country and their interaction

## Metadata for "Burdon_et_al_Fig.3_catchment_riparian_data.csv" (Numbers = column)
## 1.   Site_number: Site number at the country level
## 2.   Country: Country where data was collected (Sweden, Norway, Romania, Belgium)	
## 3.   Type_1:	Site type (forest, unbuffered, buffered)
## 4.   Type_2:	Site type (paired or reference)
## 5.   Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 6.   Catchment_PC1: Catchment impact PC1 (42.3%)	
## 7.   Catchment_PC2: Catchment impact PC2 (17.5%)	
## 8.   Riparian_PC1:	Riparian condition PC1 (36.7%)
## 9.   Riparian_PC2:	Riparian condition PC2 (14.2%)

### Load R packages for data manipulation
library(tidyr)

## Managing data frames
library(plyr)
library(dplyr)

### Load R packages for community analysis
require(vegan)

### For data visualisation
library(gplots)
require(ggplot2)
require(grid)
require(Cairo)
require(ggord)
require(ggExtra)
require(ggpubr)
require(grid)

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
library(emmeans)
library(readr)
library(stringr)

#############################################################################################
##                                                                                         ##
##                                       LOAD DATA                                         ##        
##                                                                                         ##    
#############################################################################################

# Load data
Catchment_riparian_data  <- read_csv("Burdon_et_al_Fig.3_catchment_riparian_data.csv")

#############################################################################################
##                                                                                         ##
##                          CREATE BOXPLOTS OF THE PCA GRADIENTS                           ##                
##                                                                                         ##    
#############################################################################################

## Change the order of data so that they are consistent
Catchment_riparian_data$Country <- factor(Catchment_riparian_data$Country, levels=c("SE","NO","RO","BE"))
Catchment_riparian_data$Type_1 <- factor(Catchment_riparian_data$Type_1, levels=c("UBF","FBF","FOR"))

# New facet label names for Treatment
Country.labs <- c("Sweden","Norway","Romania","Belgium")
names(Country.labs) <- c("SE","NO","RO","BE")

pCatchment <- ggplot(data=Catchment_riparian_data) +
  geom_boxplot(data=Catchment_riparian_data, aes(y=Catchment_PC1, x=Type_1, colour=Country, fill=Country), alpha=0.6) +
  theme_bw() +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country",
                      guide = FALSE) +
  # Create fill colours for the ellipses
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    name = "Country",
                    guide = FALSE
  ) +
  labs(x="Site type", y="Catchment Impact PC1") +
  theme(legend.position="left", legend.box = "vertical",
        axis.title = element_blank(),
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 14, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 14, face="bold", color="black"),
        axis.text.y=element_text(size = 11, color="black"),
        axis.text.x=element_text(size = 11, color="black", angle = 90, vjust = 0.5, hjust=1)
  ) +
        facet_wrap(. ~ Country,ncol=4,
                   labeller = labeller(Country = Country.labs))


pRiparian <- ggplot(data=Catchment_riparian_data) +
  geom_boxplot(data=Catchment_riparian_data, aes(y=Riparian_PC1, x=Type_1, colour=Country, fill=Country), alpha=0.6) +
  theme_bw() +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country",
                      guide = FALSE) +
  # Create fill colours for the ellipses
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    name = "Country",
                    guide = FALSE
  ) +
  labs(x="Site type", y="Riparian Condition PC1") +
  theme(legend.position="left", legend.box = "vertical",
        axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 14, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 14, face="bold", color="black"),
        axis.text.y=element_text(size = 11, color="black"),
        axis.text.x=element_text(size = 11, color="black", angle = 90, vjust = 0.5, hjust=1)
  ) +
  facet_wrap(. ~ Country,ncol=4,
             labeller = labeller(Country = Country.labs))


##***********************************************  
##          CREATE BOX PLOT: Fig.3  
##***********************************************

## Create multiple
png(filename="Burdon_et_al_Fig.3.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=5, 
    pointsize=20, 
    res=600)

ggarrange(pCatchment, pRiparian, ncol=2, labels = c("(a)","(b)"))

dev.off()


#############################################################################################
##                                                                                         ##
##                        TEST DIFFERENCES IN THE PCA GRADIENTS                            ##                
##                                                                                         ##    
#############################################################################################

##***********************************************  
##              CATCHMENT PC1 
##***********************************************

## Test differences in catchment gradient
M1.1 <- blmer(Catchment_PC1 ~ Type_1*Country + (1|Site_block), data = Catchment_riparian_data, REML=T,
            control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M1.1)
car::Anova(M1.1,"III",test.statistic=c("F"))

## Contrast by Site type
S.I <- lsmeans(M1.1, ~ Type_1)
pairs(S.I, adjust = "tukey")

## Contrast by Country
S.I <- lsmeans(M1.1, ~ Country)
pairs(S.I, adjust = "tukey")

## Contrast by Country and Site type
S.I <- lsmeans(M1.1, ~ Type_1 | Country)
pairs(S.I, adjust = "tukey")

## Contrast by Site type and Country
S.I <- lsmeans(M1.1, ~ Country | Type_1)
pairs(S.I, adjust = "tukey")


## Model testing Type only without Country as fixed or mixed effect
M1.2 <- blmer(Catchment_PC1 ~ Type_1 + (1|Country/Site_block), data = Catchment_riparian_data, REML=T,
              control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M1.2)
car::Anova(M1.2,"II",test.statistic=c("F"))

## Contrast by Site type
means.Type_1 <- lsmeans(M1.2, specs = "Type_1")
pairs(means.Type_1, adjust = "tukey")


## Model testing Type only without Country as fixed or mixed effect
M1.3 <- lmer(Catchment_PC1 ~ Type_1 + (1|Site_block), data = Catchment_riparian_data, REML=T,
             control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M1.3)
car::Anova(M1.3,"II",test.statistic=c("F"))

## Contrast by Site type
S.I <- lsmeans(M1.3, ~ Type_1)
pairs(S.I, adjust = "tukey")


##***********************************************  
##                RIPARIAN PC1 
##***********************************************

## Test differences in riparian PCA gradient
M2.1 <- blmer(Riparian_PC1 ~ Type_1*Country + (1|Site_block), data = Catchment_riparian_data, REML=T,
            control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M2.1)
car::Anova(M2.1,"III",test.statistic=c("F"))

## Contrast by Site type
S.I <- lsmeans(M2.1, ~ Type_1)
pairs(S.I, adjust = "tukey")

## Contrast by Country
S.I <- lsmeans(M2.1, ~ Country)
pairs(S.I, adjust = "tukey")

## Contrast by Country and Site type
S.I <- lsmeans(M2.1, ~ Type_1 | Country)
pairs(S.I, adjust = "tukey")

## Contrast by Site type and Country
S.I <- lsmeans(M2.1, ~ Country | Type_1)
pairs(S.I, adjust = "tukey")


## Model testing Type only without Country as fixed or mixed effect
M2.2 <- blmer(Riparian_PC1 ~ Type_1 + (1|Country/Site_block), data = Catchment_riparian_data, REML=T,
              control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M2.2)
car::Anova(M2.2,"II",test.statistic=c("F"))

## Contrast by Site type
means.Type_1 <- lsmeans(M2.2, specs = "Type_1")
pairs(means.Type_1, adjust = "tukey")


## Model testing Type only without Country as fixed or mixed effect
M2.3 <- lmer(Riparian_PC1 ~ Type_1 + (1|Site_block), data = Catchment_riparian_data, REML=T,
             control = lmerControl(optimizer ="Nelder_Mead"))
sjPlot::tab_model(M2.3)
car::Anova(M2.3,"II",test.statistic=c("F"))

## Contrast by Site type
S.I <- lsmeans(M2.3, ~ Type_1)
pairs(S.I, adjust = "tukey")






