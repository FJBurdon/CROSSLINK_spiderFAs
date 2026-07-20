###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Fig.8 in "Fatty acid biomarkers reveal landscape influences on                ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 4 January 2021
# Updated: 7 June 2025 

## This R script performs the structural equation modelling (SEM) for Figs.8 and S8 in Burdon et al. (in press) Ecological Monographs
## The main goal is to assess the correlation between realized trophic connectivity in riparian spiders (indicated by fatty acids) by
## accounting for key contingencies related to spider community composition, prey aquatic dispersal traits, algal productivity, and 
## environmental context (ecosystem size indicated by stream width, riparian condition PC1, and catchment impact PC1)
## This script also generates the Table 5 results reported in Burdon et al. (in press) Ecological Monographs

## Metadata for "Burdon_et_al_Fig.8_multivariate_data.csv" (Numbers = column)
## 1.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 2.     Taxa_code: Short name for arachnid family or taxonomic group	
## 3.     Type_1:	Site type (forest, unbuffered, buffered)
## 4.     Type_2:	Site type (paired or reference)
## 5.     Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 6.     Site_no: Site number at the country level
## 7.     Site_name: Site name
## 8.     Group: Arachnida	
## 9.     Genus_family: Family or taxonomic groups used for arachnids	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest level: Lowest taxonomic level used for arachnids
## 12.    Tribe
## 13-49. FA1-FA37:	Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
## 50.    EPA_perc: % concentration of Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 51.    DHA_perc: % concentration of Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 52.    LNA_perc: % concentration of Linoleic acid (LNA), an omega-6 fatty acid
## 53.    ALA_perc: % concentration of Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 54.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 55.    chl_a_m2_day: Algal accrual on tiles
## 56.    Width: Stream wetted channel width (m)
## 57.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 58.    Catchment_PC1: Catchment impact PC1 (42.3%)	
## 59.     s3: Maximum body size (cm) ≥ 0.5–1
## 60.    cd1: Life cycle duration ≤ 1 year
## 61.    cy2: Potential number of reproductive cycles per year - One
## 62.    dis3: Dispersal strategy Aerial passive
## 63.    dis4: Dispersal strategy Aerial active
## 64.    life1: Adult life span < 1 week
## 65.    life2: Adult life span ≥ 1 week – 1 month
## 66.    wnb2: Wing pair type - 1 pair + halters
## 67.    wnb3: Wing pair type - 1 pair + 1 pair of small hind wings
## 68.    wnb5: Wing pair type - 2 similar-sized pairs
## 69.    EPA_ALA_ratio: Ratio of EPA to ALA

## Metadata for "Burdon_et_al_Fig.8_Lycosidae_SEM_data_std" (Numbers = column)
## Same as "Burdon_et_al_Fig.8_multivariate_data.csv" except Tribe column has been removed
## And a column has been added for Lycosidae catch per unit effort (individuals per m2 per min) 
## NOTE: data is already transformed and standardised

## Metadata for "Burdon_et_al_Fig.8_Tetragnathidae_SEM_data_std.csv" (Numbers = column)
## Same as "Burdon_et_al_Fig.8_multivariate_data.csv" except Tribe column has been removed
## And a column has been added for Tetragnathidae catch per unit effort (individuals per m2 per min) 
## NOTE: data is already transformed and standardised

### Load R packages for data management
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(hms)

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
library(effects)
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
library(optimx)

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

## Load SEM data - includes PCA gradients etc.
Spider_SEM_data <- read_csv("Burdon_et_al_Fig.8_multivariate_data.csv")

## Load SEM data for Spider Families - already transformed and standardized
Spider_SEM_data_Lycos <- read_csv("Burdon_et_al_Fig.8_Lycosidae_SEM_data_std.csv")
Spider_SEM_data_Tetra <- read_csv("Burdon_et_al_Fig.8_Tetragnathidae_SEM_data_std.csv")

#############################################################################################
##                                                                                         ##
##                                    SUBSET DATA INTO GROUPS                              ##        
##                                                                                         ##    
#############################################################################################

## Create subset data for further analysis

## Spider Mode of Hunting
Spider_SEM_data_ground <- Spider_SEM_data[-which(Spider_SEM_data$Mode=="Web"),]  
Spider_SEM_data_web <- Spider_SEM_data[-which(Spider_SEM_data$Mode=="Ground"),]

#############################################################################################
##                                                                                         ##
##                                   SEM - ALL SPIDERS                                     ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                              TRANSFORM AND STANDARDIZE DATA                             ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset for transformations
Spider_SEM_data_tsfd  <- Spider_SEM_data

## Transform percentages using logit transformation
Spider_SEM_data_tsfd[,c(13:53)] <- logit(Spider_SEM_data_tsfd[,c(13:53)]/100)

## Check selection
colnames(Spider_SEM_data_tsfd)

### PUFAs: Data transformations before analysis
Spider_SEM_data_tsfd [,c(56)] <- log(Spider_SEM_data_tsfd [,c(56)])#log transform Width
Spider_SEM_data_tsfd [,c(55)] <- log1p(Spider_SEM_data_tsfd [,c(55)])#log + 1 transform Chl-a
Spider_SEM_data_tsfd [,c(54,59:68)] <- log1p(Spider_SEM_data_tsfd [,c(54,59:68)])#log + 1 transform traits

## Remove tribe since not applicable to arachnid data
Spider_SEM_data_std <- Spider_SEM_data_tsfd[,-12]

## Check column names
colnames(Spider_SEM_data_std[,c(49,50:68)])

## Standardize data
Spider_SEM_data_std[,c(49,50:68)] <- decostand(Spider_SEM_data_std[,c(49,50:68)],"standardize",na.action=na.exclude)

## Remove incomplete records
Spider_SEM_data_std <- na.omit(Spider_SEM_data_std)

#############################################################################################
##                                                                                         ##
##                                  SEM FOR EPA:ALA RATIO                                  ## 
##                                     DIS4 AND DIS3                                       ##
##                                                                                         ##    
#############################################################################################

# Create component models and store in list
PUFA_pSEM_F1 = psem(
  
  # Predicting connectivity response
  Spider_response = lmer(EPA_ALA_ratio ~ dis4 + dis3 + Width + Riparian_PC1 + (1|Site_block/Site_name) + (1|Genus_family),
                         control = lmerControl(optimizer = "Nelder_Mead"),
                         data = Spider_SEM_data_std),
  
  # Predicting spider response
  Spider_response = lmer(Spider_body_size ~ Riparian_PC1 + Width + (1|Site_block),
                         control = lmerControl(optimizer = "Nelder_Mead"),
                         data = Spider_SEM_data_std),
  
  # Predicting stream invert traits response
  Insect_response1 = lmer(dis4 ~ chl_a_m2_day + Riparian_PC1 + Catchment_PC1 + (1|Site_block),
                          control = lmerControl(optimizer = "Nelder_Mead"),
                          data = Spider_SEM_data_std),
  
  # Predicting stream invert traits response
  Insect_response2 = lmer(dis3 ~ Catchment_PC1 + (1|Site_block),
                          control = lmerControl(optimizer = "Nelder_Mead"),
                          data = Spider_SEM_data_std),
  
  # Predicting basal resource response
  Algae_response = lmer(chl_a_m2_day ~ dis3 + Riparian_PC1 + (1|Country/Site_block),
                        control = lmerControl(optimizer = "Nelder_Mead"),
                        data = Spider_SEM_data_std))

# Adding correlated error terms
PUFA_pSEM_F1  <- update(PUFA_pSEM_F1,
                        dis4 %~~% dis3,
                        dis4 %~~% Spider_body_size,
                        EPA_ALA_ratio %~~% Spider_body_size
)

## Information criterion - AICc
AIC(PUFA_pSEM_F1, aicc = T)

# Run goodness-of-fit tests
summary(PUFA_pSEM_F1, standardize = "scale", intercepts = FALSE)

## AIC + AICc Final model
## 2271.025 2272.973 30 238

## AIC + AICc Alternative model
## EPA_ALA_ratio ~  chl_a_m2_day...
## 2275.418 2277.545 31 238

## AIC + AICc dis4 only
## 2267.475 2269.261 29 238

## AIC + AICc dis3 only
## 2275.311 2277.097 29 238

# Explore individual model fits
rsquared(PUFA_pSEM_F1)

# Evaluate path significance using standardized coefficients
coefs(PUFA_pSEM_F1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(PUFA_pSEM_F1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_PUFA_F1_dis4_dis3_SpiderAll_EPA_ALA_ratio_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(PUFA_pSEM_F1)
Goodness_of_fit <- rsquared(PUFA_pSEM_F1)
write.csv(Goodness_of_fit, "SEM_PUFA_F1_dis4_dis3_SpiderAll_EPA_ALA_ratio_rsquare.csv",row.names = T)

#############################################################################################
##                                                                                         ##
##                                    KEY SPIDER FAMILIES                                  ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                                   SEM - TETRAGNATHIDAE                                  ##        
##                                                                                         ##    
#############################################################################################

## Remove incomplete records
Spider_SEM_data_Tetra <- na.omit(Spider_SEM_data_Tetra)

#############################################################################################
##                                                                                         ##
##                                  SEM FOR EPA:ALA RATIO                                  ## 
##                                      DIS4 AND DIS3                                      ##
##                                                                                         ##    
#############################################################################################

## Compare the model used for EPA (dis4 + dis3)
## Create component models and store in list
PUFA_pSEM_T1 = psem(
  
  # Predicting connectivity response
  Spider_response = lm(EPA_ALA_ratio ~ dis4 + dis3 + Width,
                       data = Spider_SEM_data_Tetra),
  
  # Predicting spider response
  Spider_response = lm(Tetragnathidae ~ 1, 
                       data = Spider_SEM_data_Tetra),
  
  # Predicting stream invert traits response
  Insect_response2 = lm(dis4 ~ Catchment_PC1 + chl_a_m2_day,
                        data = Spider_SEM_data_Tetra),
  
  # Predicting stream invert traits response
  Insect_response3 = lm(dis3 ~ Catchment_PC1 + Riparian_PC1,
                        data = Spider_SEM_data_Tetra),
  
  # Predicting basal resource response
  Algae_response = lmer(chl_a_m2_day ~ dis3 + Riparian_PC1 + (1|Country),  
                        data = Spider_SEM_data_Tetra,
                        REML = T, 
                        control = lmerControl(
                          optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))))

PUFA_pSEM_T1  <- update(PUFA_pSEM_T1,
                        dis4 %~~% dis3)

## Information criterion - AIC and AICc
AIC(PUFA_pSEM_T1, aicc = T)
## 708.804 713.263 20 53

# Run goodness-of-fit tests
summary(PUFA_pSEM_T1, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(PUFA_pSEM_T1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(PUFA_pSEM_T1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_PUFA_T1_dis4_dis3_Tetragnathid_EPA_ALA_ratio_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(PUFA_pSEM_T1)
Goodness_of_fit <- rsquared(PUFA_pSEM_T1)
write.csv(Goodness_of_fit, "SEM_PUFA_T1_dis4_dis3_Tetragnathid_EPA_ALA_ratio_rsquare.csv",row.names = T)


#############################################################################################
##                                                                                         ##
##                                   SEM - LYCOSIDAE                                       ##        
##                                                                                         ##    
#############################################################################################

## Remove incomplete records
Spider_SEM_data_Lycos <- na.omit(Spider_SEM_data_Lycos)

#############################################################################################
##                                                                                         ##
##                                  SEM FOR EPA:ALA RATIO                                  ## 
##                                      DIS4 AND DIS3                                      ##
##                                                                                         ##    
#############################################################################################

## Compare the model used for EPA (dis4 + dis3)
## Create component models and store in list
PUFA_pSEM_L1 = psem(
  
  # Predicting connectivity response
  Spider_response = lm(EPA_ALA_ratio ~ dis4 + dis3 + Width + Catchment_PC1,
                       data = Spider_SEM_data_Lycos),
  
  # Predicting spider response
  Spider_response = lm(Lycosidae ~ Riparian_PC1 + dis4, 
                       data = Spider_SEM_data_Lycos),
  
  # Predicting stream invert traits response
  Insect_response2 = lm(dis4 ~ Catchment_PC1 + chl_a_m2_day,
                        data = Spider_SEM_data_Lycos),
  
  # Predicting stream invert traits response
  Insect_response3 = lm(dis3 ~ Catchment_PC1 + Riparian_PC1,
                        data = Spider_SEM_data_Lycos),
  
  # Predicting basal resource response
  Algae_response = lmer(chl_a_m2_day ~ dis3 + (1|Country),  
                        data = Spider_SEM_data_Lycos,
                        REML = T, 
                        control = lmerControl(
                          optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))))

PUFA_pSEM_L1  <- update(PUFA_pSEM_L1,
                        dis4 %~~% dis3)

## Information criterion - AIC and AICc
AIC(PUFA_pSEM_L1, aicc = T)
## 819.994 824.328 22 6

# Run goodness-of-fit tests
summary(PUFA_pSEM_L1, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(PUFA_pSEM_L1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(PUFA_pSEM_L1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_PUFA_L1_dis4_dis3_Lycosid_EPA_ALA_ratio_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(PUFA_pSEM_L1)
Goodness_of_fit <- rsquared(PUFA_pSEM_L1)
write.csv(Goodness_of_fit, "SEM_PUFA_L1_dis4_dis3_Lycosid_EPA_ALA_ratio_rsquare.csv",row.names = T)

#############################################################################################
##                                                                                         ##
##                                SEM - WEB-BUILDING SPIDERS                               ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                              TRANSFORM AND STANDARDIZE DATA                             ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset for transformations
Spider_SEM_data_web_tsfd  <- Spider_SEM_data_web

## Transform percentages using logit transformation
Spider_SEM_data_web_tsfd[,c(13:53)] <- logit(Spider_SEM_data_web_tsfd[,c(13:53)]/100)

## Check selection
colnames(Spider_SEM_data_web_tsfd)

### PUFAs: Data transformations before analysis
Spider_SEM_data_web_tsfd [,c(56)] <- log(Spider_SEM_data_web_tsfd [,c(56)])#log transform Width
Spider_SEM_data_web_tsfd [,c(55)] <- log1p(Spider_SEM_data_web_tsfd [,c(55)])#log + 1 transform Chl-a
Spider_SEM_data_web_tsfd [,c(54,59:68)] <- log1p(Spider_SEM_data_web_tsfd [,c(54,59:68)])#log + 1 transform traits

## Remove tribe since not applicable to arachnid data
Spider_SEM_data_web_std <- Spider_SEM_data_web_tsfd[,-12]

## Check column names
colnames(Spider_SEM_data_web_std[,c(49,50:68)])

## Standardize data
Spider_SEM_data_web_std[,c(49,50:68)] <- decostand(Spider_SEM_data_web_std[,c(49,50:68)],"standardize",na.action=na.exclude)

## Remove incomplete records
Spider_SEM_data_web_std <- na.omit(Spider_SEM_data_web_std)

#############################################################################################
##                                                                                         ##
##                                  SEM FOR EPA:ALA RATIO                                  ## 
##                                     DIS4 AND DIS3                                       ##
##                                                                                         ##    
#############################################################################################

# Create component models and store in list
PUFA_pSEM_W1 = psem(
  
  # Predicting connectivity response
  Spider_response = lmer(EPA_ALA_ratio ~ dis4 + dis3 + Width + (1|Genus_family),
                         control = lmerControl(optimizer = "Nelder_Mead"),
                         data = Spider_SEM_data_web_std),
  
  # Predicting spider response
  Spider_response = lm(Spider_body_size ~ Riparian_PC1,
                       data = Spider_SEM_data_web_std),
  
  # Predicting stream invert traits response
  Insect_response1 = lm(dis4 ~ Catchment_PC1 + chl_a_m2_day,
                        data = Spider_SEM_data_web_std),
  
  # Predicting stream invert traits response
  Insect_response2 = lm(dis3 ~ Catchment_PC1,
                        data = Spider_SEM_data_web_std),
  
  # Predicting basal resource response
  Algae_response = lmer(chl_a_m2_day ~ dis3 + Riparian_PC1 + (1|Country/Site_block),
                        control = lmerControl(optimizer = "Nelder_Mead"),
                        data = Spider_SEM_data_web_std))

PUFA_pSEM_W1  <- update(PUFA_pSEM_W1,
                        dis4 %~~% dis3,
                        EPA_ALA_ratio %~~% Spider_body_size
)

## Information criterion - AIC and AICc
AIC(PUFA_pSEM_W1, aicc = T)
## 1635.607 1637.609 22 134

# Run goodness-of-fit tests
summary(PUFA_pSEM_W1, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(PUFA_pSEM_W1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(PUFA_pSEM_W1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_PUFA_W1_dis4_dis3_Web_spider_EPA_ALA_ratio_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(PUFA_pSEM_W1)
Goodness_of_fit <- rsquared(PUFA_pSEM_W1)
write.csv(Goodness_of_fit, "SEM_PUFA_W1_dis4_dis3_Web_spider_EPA_ALA_ratio_rsquare.csv",row.names = T)

#############################################################################################
##                                                                                         ##
##                               SEM - GROUND-HUNTING SPIDERS                              ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                              TRANSFORM AND STANDARDIZE DATA                             ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset for transformations
Spider_SEM_data_ground_tsfd  <- Spider_SEM_data_ground

## Transform percentages using logit transformation
Spider_SEM_data_ground_tsfd[,c(13:53)] <- logit(Spider_SEM_data_ground_tsfd[,c(13:53)]/100)

## Check selection
colnames(Spider_SEM_data_ground_tsfd)

### PUFAs: Data transformations before analysis
Spider_SEM_data_ground_tsfd [,c(56)] <- log(Spider_SEM_data_ground_tsfd [,c(56)])#log transform Width
Spider_SEM_data_ground_tsfd [,c(55)] <- log1p(Spider_SEM_data_ground_tsfd [,c(55)])#log + 1 transform Chl-a
Spider_SEM_data_ground_tsfd [,c(54,59:68)] <- log1p(Spider_SEM_data_ground_tsfd [,c(54,59:68)])#log + 1 transform traits

## Remove tribe since not applicable to arachnid data
Spider_SEM_data_ground_std <- Spider_SEM_data_ground_tsfd[,-12]

## Check column names
colnames(Spider_SEM_data_ground_std[,c(49,50:68)])

## Standardize data
Spider_SEM_data_ground_std[,c(49,50:68)] <- decostand(Spider_SEM_data_ground_std[,c(49,50:68)],"standardize",na.action=na.exclude)

## Remove incomplete records
Spider_SEM_data_ground_std <- na.omit(Spider_SEM_data_ground_std)

#############################################################################################
##                                                                                         ##
##                                  SEM FOR EPA:ALA RATIO                                  ## 
##                                     DIS4 AND DIS3                                       ##
##                                                                                         ##    
#############################################################################################

# Create component models and store in list
PUFA_pSEM_G1 = psem(
  
  # Predicting connectivity response
  Spider_response = lmer(EPA_ALA_ratio ~ dis4 + dis3 + Width + Catchment_PC1 + (1|Genus_family),
                         control = lmerControl(optimizer = "Nelder_Mead"),
                         data = Spider_SEM_data_ground_std),
  
  # Predicting spider response
  Spider_response = lm(Spider_body_size ~ Riparian_PC1,
                       data = Spider_SEM_data_ground_std),
  
  # Predicting stream invert traits response
  Insect_response1 = lm(dis4 ~ Catchment_PC1 + chl_a_m2_day,
                        data = Spider_SEM_data_ground_std),
  
  # Predicting stream invert traits response
  Insect_response2 = lm(dis3 ~ Catchment_PC1 + Riparian_PC1,
                        data = Spider_SEM_data_ground_std),
  
  # Predicting basal resource response
  Algae_response = lmer(chl_a_m2_day ~ dis3 + (1|Country/Site_block),
                        control = lmerControl(optimizer = "Nelder_Mead"),
                        data = Spider_SEM_data_ground_std))

PUFA_pSEM_G1  <- update(PUFA_pSEM_G1,
                        dis4 %~~% dis3,
                        EPA_ALA_ratio %~~% Spider_body_size
)

## Information criterion - AIC and AICc
AIC(PUFA_pSEM_G1, aicc = T)
## 1266.664 1269.491 23 104

# Run goodness-of-fit tests
summary(PUFA_pSEM_G1, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(PUFA_pSEM_G1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(PUFA_pSEM_G1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_PUFA_G1_dis4_dis3_Ground_spider_EPA_ALA_ratio_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(PUFA_pSEM_G1)
Goodness_of_fit <- rsquared(PUFA_pSEM_G1)
write.csv(Goodness_of_fit, "SEM_PUFA_G1_dis4_dis3_Ground_spider_EPA_ALA_ratio_rsquare.csv",row.names = T)


