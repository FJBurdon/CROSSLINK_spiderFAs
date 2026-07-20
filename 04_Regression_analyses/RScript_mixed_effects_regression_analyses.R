###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Fig.7 in "Fatty acid biomarkers reveal landscape influences on                ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 4 July 2021
# Updated: 7 June 2025 

## This R script performs the linear mixed effects modelling for Figs.7, S5, and S6 in Burdon et al. (in press) Ecological Monographs
## The main goal is to assess the correlation between realized trophic connectivity in riparian spiders (indicated by fatty acids) and
## community-weighted mean trait abundances of stream invertebrates including aerial active and aerial passive dispersal modalities
## This script also generates the results reported in Table 4 and S9 in Burdon et al. (in press) Ecological Monographs

## Metadata for "Burdon_et_al_Fig.7_multivariate_data.csv" (Numbers = column)
## 1.     Site_no: Site number	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Taxa_code: Short name for arachnid family or taxonomic group	
## 4.     Genus_family: Family or taxonomic groups used for arachnids	
## 5.     Site_name: Site name
## 6.     Type_1:	Site type (forest, unbuffered, buffered)
## 7.     Type_2:	Site type (paired or reference)
## 8.     Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 9.     Group: Arachnida	
## 10.    Mode: Mode of hunting for arachnids
## 11.    EPA: concentration of EPA mg/g, Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 12.    Ratio_EPA_LNA: Ratio of EPA to Linoleic acid (LNA), an omega-6 fatty acid
## 13.    EPA_perc: %EPA
## 14.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 15.     s3: Maximum body size (cm) ≥ 0.5–1
## 16.    cd1: Life cycle duration ≤ 1 year
## 17.    dis4: Dispersal strategy Aerial active
## 18.    wnb3: Wing pair type - 1 pair + 1 pair of small hind wings
## 19.    Insect_perc: proportion of invertebrates in Surber sample that are insects
## 20.    Insect_m2: insect density per square-meter from Surber sample
## 21.    EPT_perc: proportion of invertebrates in Surber sample that are EPT insects 
## 22.    EPT_m2: EPT insect density per square-meter from Surber sample
## 23.    Chironomid_perc: proportion of invertebrates in Surber sample that are Chironomidae
## 24.    Chironomid_m2: Chironomidae density per square-meter from Surber sample
## 25.    chl_a_m2_day: Algal accrual on tiles
## 26.    Width: Stream wetted channel width (m)
## 27.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 28.    Riparian_RCI: Riparian Condition Index, see Burdon et al. 2020 https://doi.org/10.3390/w12041178
## 29.    Catchment_PC1: Catchment impact PC1 (42.3%)	
## 30.    DHA: concentration mg/g of DHA, Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 31.    LNA: concentration mg/g of LNA, Linoleic acid (LNA), an omega-6 fatty acid
## 32.    ALA: concentration mg/g of ALA, Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 33.    DHA_perc: %DHA  
## 34.    LNA_perc: %LNA
## 35.    ALA_perc: #ALA
## 36.    EPA_ALA_ratio: Ratio of EPA to ALA
## 37.    dis3: Dispersal strategy Aerial passive
## 38.    Body_Size:	Mean body size of spider Family sampled for FAs

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
Spider_LME_data <- read_csv("Burdon_et_al_Fig.7_multivariate_data.csv")
 
#############################################################################################
##                                                                                         ##
##                                     TRANSFORMATIONS                                     ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset for transformations
Spider_LME_data_tsfd <- Spider_LME_data

## LOG
#colnames(Spider_LME_data_tsfd[,c(11,20,26,28,38)])
# "EPA"          "Insect_m2"    "Width"        "Riparian_RCI" "Body_Size"   
Spider_LME_data_tsfd[,c(11,20,26,28,38)] <- log(Spider_LME_data_tsfd[,c(11,20,26,28,38)])#log transform PUFA

## LOG+1 TRAITS
#  colnames(Spider_LME_data_tsfd[,c(14:18,37)])
#  "Spider_body_size"   "s3"   "cd1"   "dis4"    "wnb3"  "dis3"     
Spider_LME_data_tsfd[,c(12,14:18,37)] <- log1p(Spider_LME_data_tsfd[,c(12,14:18,37)])#sqrt transform traits

## LOGIT
#colnames(Spider_LME_data_tsfd[,c(13,19,21,23)])
# "EPA_perc"        "Insect_perc"     "EPT_perc"        "Chironomid_perc"
Spider_LME_data_tsfd[,c(13,19,21,23)] <- logit(Spider_LME_data_tsfd[,c(13,19,21,23)]/100)#logit % Arable cropping

## LOG+1
#colnames(Spider_LME_data_tsfd[,c(22,24,25)])
# "EPT_m2"    "Chironomid_m2"    "chl_a_m2_day" 
Spider_LME_data_tsfd[,c(22,24,25)] <- log1p(Spider_LME_data_tsfd[,c(22,24,25)])#log + 1 transform
## Transform additional data added to dataframe

## LOG
#colnames(Spider_LME_data_tsfd[,c(30:32,36)])
# "DHA"           "LNA"           "ALA"    "EPA_ALA_ratio"
Spider_LME_data_tsfd[,c(30:32,36)] <- log(Spider_LME_data_tsfd[,c(30:32,36)])#log transform PUFA concentrations

## ASIN SQRT
#colnames(Spider_LME_data_tsfd[,c(33:35)])
# "DHA_perc"      "LNA_perc"      "ALA_perc"
Spider_LME_data_tsfd[,c(33:35)] <- asin(sqrt(Spider_LME_data_tsfd[,c(33:35)]/100))

#############################################################################################
##                                                                                         ##
##                                    SUBSET DATA INTO GROUPS                              ##        
##                                                                                         ##    
#############################################################################################

## Create subset data for further analysis
## Spider Mode of Hunting
Spider_LME_data_ground_tsfd <- Spider_LME_data_tsfd[-which(Spider_LME_data_tsfd$Mode=="Web"),]  
Spider_LME_data_web_tsfd <- Spider_LME_data_tsfd[-which(Spider_LME_data_tsfd$Mode=="Ground"),]

## Spider Families
Spider_LME_data_Lycos_tsfd <- Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Genus_family=="Lycosidae"),]  
Spider_LME_data_Tetra_tsfd <- Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Genus_family=="Tetragnathidae"),]  

#############################################################################################
##                                                                                         ##
##                                    STANDARDIZATION                                      ##        
##                                                                                         ##    
#############################################################################################

## Now remove all missing values Rows 250 -> 234 (i.e., 16 rows with missing values)
Spider_LME_data_std <- na.omit(Spider_LME_data_tsfd)
Spider_LME_data_ground_std <- na.omit(Spider_LME_data_ground_tsfd)
Spider_LME_data_web_std <- na.omit(Spider_LME_data_web_tsfd)
Spider_LME_data_Lycos_std <- na.omit(Spider_LME_data_Lycos_tsfd)
Spider_LME_data_Tetra_std <- na.omit(Spider_LME_data_Tetra_tsfd)

## Standardise data for further analysis
Spider_LME_data_std[,c(11:38)] <- decostand(Spider_LME_data_std[,c(11:38)],"standardize",na.rm = TRUE)
Spider_LME_data_web_std[,c(11:38)] <- decostand(Spider_LME_data_web_std[,c(11:38)],"standardize",na.rm = TRUE)
Spider_LME_data_ground_std[,c(11:38)] <- decostand(Spider_LME_data_ground_std[,c(11:38)],"standardize",na.rm = TRUE)
Spider_LME_data_Lycos_std[,c(11:38)] <- decostand(Spider_LME_data_Lycos_std[,c(11:38)],"standardize",na.rm = TRUE)
Spider_LME_data_Tetra_std[,c(11:38)] <- decostand(Spider_LME_data_Tetra_std[,c(11:38)],"standardize",na.rm = TRUE)

## Relabel row numbers sequentially
row.names(Spider_LME_data_std) <- NULL

## Check site headers
head(Spider_LME_data_std)

## Relabel the Forest sites 
Spider_LME_data_std[which(Spider_LME_data_std$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_std[which(Spider_LME_data_std$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_std[which(Spider_LME_data_std$Type_1=="LDS"),6]   <- "FOR"

#############################################################################################
##                                                                                         ##
##                         LME Correlation EPA% - ALL SPIDERS DIS4                         ##        
##                                                                                         ##    
#############################################################################################

## Relabel the Forest sites 
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="LDS"),6]   <- "FOR"
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_2=="LDS"),7]   <- "FOR"

## Check this worked
unique(Spider_LME_data_tsfd$Type_1)
unique(Spider_LME_data_tsfd$Type_2)

## Create model with transformed data for plot
M1 <- blmer(EPA_perc ~ dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_tsfd,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_dis4 <- effects::effect(term= "dis4", mod = M1)
summary(effects_dis4) #output of what the values are

# Save the effects values as a df:
x_dis4 <- as.data.frame(effects_dis4)

## Check unique
unique(Spider_LME_data_tsfd$Type_1)

## Head
head(Spider_LME_data_tsfd)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Create plot
png(filename="Burdon_et_al_Fig.7.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=20, 
    res=600)

dis4_plot <- ggplot() + 
  #1
  #stat_smooth(data=Spider_LME_data_tsfd, aes(dis4, EPA_perc), method = "lm", color="black", se=T) +
  geom_line(data=x_dis4, aes(x=dis4, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_dis4, aes(x=dis4, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Spider_LME_data_tsfd, aes(dis4, EPA_perc, colour=Country, shape=Type_1), size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  scale_shape_manual(breaks=c("UBF", "FBF", "FOR"),
                     values=c(15, 17, 16),
                     labels = c("Unbuffered", "Buffered", "Forest"),
                     name = "Site type"
                     #, guide = FALSE
  ) +
  #5
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country"
                      #, guide = FALSE
  ) +
  #6
  xlab(expression(bold("Aquatic invert. aerial active dispersal (log["*italic(x)*"+1] CWM)"))) +
  ylab("Trophic connectivity (logit spider %EPA)") +
  #
  annotate("text", x=0.4, y=-2.25, label=expression(italic(R)[m]^2*~"= 16.7%,  "*italic(R)[c]^2*~"= 59.7%"), size=6, fontface=2) + 
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title.y=element_text(vjust=1.25, size = 18, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 18, face="bold", color="black"),
        axis.text.y=element_text(size = 14, color="black"),
        axis.text.x=element_text(size = 14, color="black"))

## Create marginal plot
ggMarginal(dis4_plot,groupColour = TRUE, groupFill = TRUE)

dev.off()

#############################################################################################
##                                                                                         ##
##                         LME Correlation EPA:ALA - ALL SPIDERS DIS4                      ##        
##                                                                                         ##    
#############################################################################################


## Relabel the Forest sites 
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_tsfd[which(Spider_LME_data_tsfd$Type_1=="LDS"),6]   <- "FOR"

## Create model with transformed data for plot
M1 <- blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_tsfd,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_dis4 <- effects::effect(term= "dis4", mod = M1)
summary(effects_dis4) #output of what the values are

# Save the effects values as a df:
x_dis4 <- as.data.frame(effects_dis4)

## Check unique
unique(Spider_LME_data_tsfd$Type_1)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Create plot
png(filename="Burdon_et_al_Fig.S5.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=20, 
    res=600)

dis4_plot <- ggplot() + 
  #1
  #stat_smooth(data=Spider_LME_data_tsfd, aes(dis4, EPA_perc), method = "lm", color="black", se=T) +
  geom_line(data=x_dis4, aes(x=dis4, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_dis4, aes(x=dis4, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Spider_LME_data_tsfd, aes(dis4, EPA_ALA_ratio, colour=Country, shape=Type_1), size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  scale_shape_manual(breaks=c("UBF", "FBF", "FOR"),
                     values=c(15, 17, 16),
                     labels = c("Unbuffered", "Buffered", "Forest"),
                     name = "Site type"
                     #, guide = FALSE
  ) +
  #5
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country"
                      #, guide = FALSE
  ) +
  #6
  #labs(x="Aquatic invertebrate dispersal (CWM aerial active)", y="Trophic connectivity (Spider EPA:ALA)") +
  ylab("Trophic connectivity (log spider EPA:ALA)") +
  xlab(expression(bold("Aquatic invert. aerial active dispersal (log["*italic(x)*"+1] CWM)"))) +
  #
  annotate("text", x=0.4, y=3.2, label=expression(italic(R)[m]^2*~"= 14.9%,  "*italic(R)[c]^2*~"= 58.6%"), size=6, fontface=2) + 
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title.y=element_text(vjust=1.25, size = 18, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 18, face="bold", color="black"),
        axis.text.y=element_text(size = 14, color="black"),
        axis.text.x=element_text(size = 14, color="black"))

## Create marginal plot
ggMarginal(dis4_plot,groupColour = TRUE, groupFill = TRUE)

dev.off()

#############################################################################################
##                                                                                         ##
##                         LME Correlation EPA% - ALL SPIDERS DIS3                         ##        
##                                                                                         ##    
#############################################################################################

## Create model with transformed data for plot
M1 <- blmer(EPA_perc ~ dis3 + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_tsfd,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_dis3 <- effects::effect(term= "dis3", mod = M1)
summary(effects_dis3) #output of what the values are

# Save the effects values as a df:
x_dis3 <- as.data.frame(effects_dis3)

## Check unique
unique(Spider_LME_data_tsfd$Type_1)

## Head
head(Spider_LME_data_tsfd)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

head(Spider_LME_data_tsfd)

## Create plot
png(filename="Burdon_et_al_Fig.S6.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=20, 
    res=600)

dis3_plot <- ggplot() + 
  #1
  #stat_smooth(data=Spider_LME_data_tsfd, aes(dis4, EPA_perc), method = "lm", color="black", se=T) +
  geom_line(data=x_dis3, aes(x=dis3, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_dis3, aes(x=dis3, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Spider_LME_data_tsfd, aes(dis3, EPA_perc, colour=Country, shape=Type_1), size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  scale_shape_manual(breaks=c("UBF", "FBF", "FOR"),
                     values=c(15, 17, 16),
                     labels = c("Unbuffered", "Buffered", "Forest"),
                     name = "Site type"
                     #, guide = FALSE
  ) +
  #5
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country"
                      #, guide = FALSE
  ) +
  #6
  xlab(expression(bold("Aquatic invert. aerial passive dispersal (log["*italic(x)*"+1] CWM)"))) +
  ylab("Trophic connectivity (logit spider %EPA)") +
  #
  annotate("text", x=0.4, y=-2.25, label=expression(italic(R)[m]^2*~"= 6.7%,  "*italic(R)[c]^2*~"= 58.8%"), size=6, fontface=2) + 
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        axis.title.y=element_text(vjust=1.25, size = 18, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 18, face="bold", color="black"),
        axis.text.y=element_text(size = 14, color="black"),
        axis.text.x=element_text(size = 14, color="black"))

## Create marginal plot
ggMarginal(dis3_plot,groupColour = TRUE, groupFill = TRUE)

dev.off()

#############################################################################################
##                                                                                         ##
##                        LME Correlation EPA:ALA - ALL SPIDERS DIS3                       ##        
##                                                                                         ##    
#############################################################################################

## Create model with transformed data for plot
M1 <- blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_tsfd,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_dis3 <- effects::effect(term= "dis3", mod = M1)
summary(effects_dis3) #output of what the values are

# Save the effects values as a df:
x_dis3 <- as.data.frame(effects_dis3)

## Check unique
unique(Spider_LME_data_tsfd$Type_1)

## Head
head(Spider_LME_data_tsfd)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

head(Spider_LME_data_tsfd)


dis3_plot <- ggplot() + 
  #1
  #stat_smooth(data=Spider_LME_data_tsfd, aes(dis4, EPA_perc), method = "lm", color="black", se=T) +
  geom_line(data=x_dis3, aes(x=dis3, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_dis3, aes(x=dis3, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Spider_LME_data_tsfd, aes(dis3, EPA_ALA_ratio, colour=Country, shape=Type_1), size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  scale_shape_manual(breaks=c("UBF", "FBF", "FOR"),
                     values=c(15, 17, 16),
                     labels = c("Unbuffered", "Buffered", "Forest"),
                     name = "Site type"
                     #, guide = FALSE
  ) +
  #5
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country"
                      #, guide = FALSE
  ) +
  #6
  labs(x="Aquatic invertebrate dispersal (CWM aerial passive)", y="Trophic connectivity (Spider EPA:ALA)") +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold", color="black"),
        axis.text.y=element_text(size = 13, color="black"),
        axis.text.x=element_text(size = 13, color="black"))

## Create marginal plot
ggMarginal(dis3_plot,groupColour = TRUE, groupFill = TRUE)


#############################################################################################
##                                                                                         ##
##                          CHECK PREDICTORS - DIS4 and COUNTRY                            ##        
##                                                                                         ##    
#############################################################################################

## Take site means for Dis4
des4_mu <- Spider_LME_data_std %>% group_by(Country,Site_name,Type_1,Type_2,Site_block) %>% 
                                      summarise(dis4=mean(dis4))

## Test the correlation of dis3 and dis4
M1 <- blmer(dis4 ~ Country + (1|Site_block), 
            data=des4_mu,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1)
car::Anova(M1,"II")

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
L.S <- pairs(lsmeans(M1, ~ Country))
test(L.S, adjust = "tukey")

## Take site means for Dis3
des3_mu <- Spider_LME_data_std %>% group_by(Country,Site_name,Type_1,Type_2,Site_block) %>% 
  summarise(dis3=mean(dis3))

## Test the correlation of dis3
M2 <- blmer(dis3 ~ Country + (1|Site_block), 
            data=des3_mu,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M2)
car::Anova(M2,"II")

## Write post-hoc test output ANOVA - drop COUNTRY and nested term SITE_BLOCK
L.S <- pairs(lsmeans(M2, ~ Country))
test(L.S, adjust = "tukey")

#############################################################################################
##                                                                                         ##
##                          ALA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(ALA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(ALA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(ALA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(ALA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis3 and Country > 2

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(ALA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## See influence of interation betwwn SITE and COUNTRY
M1.6 <- blmer(ALA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## See difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##          MODE          ##
############################

## Test if Mode makes a difference
M1.7 <- blmer(ALA_perc ~ Mode*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.7) ## See correlation of fixed effects
sjPlot::tab_model(M1.7)
BIC(M1.7)
car::Anova(M1.7,"III")
performance::check_collinearity(M1.7) ## VIF > 4

## Plot data
boxplot(ALA_perc ~ Mode*Country, data=Spider_LME_data_std)
#dev.off()

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Mode | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Country | Mode))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(ALA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(ALA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## See influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),],main="Belgium")
abline(lm(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]))
## Norway
plot(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),],main="Norway")
abline(lm(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]))
## Sweden
plot(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),],main="Sweden")
abline(lm(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]))
## Romania
plot(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),],main="Romania")
abline(lm(ALA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.7,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country < 2

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##          MODE          ##
############################

## Test if Mode makes a difference
M1.7 <- blmer(EPA_perc ~ Mode*Country + (1|Site_block/Site_name) + (1|Genus_family), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
summary(M1.7) ## See correlation of fixed effects
sjPlot::tab_model(M1.7)
BIC(M1.7)
car::Anova(M1.7,"III")
performance::check_collinearity(M1.7) ## VIF > 4

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Mode | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Country | Mode))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(EPA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
            data=Spider_LME_data_std,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## See influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]))
## Norway
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),],main="Norway")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]))
## Sweden
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]))
## Romania
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),],main="Romania")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.7,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          DHA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(DHA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(DHA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(DHA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(DHA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(DHA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##          MODE          ##
############################

## Test if Mode makes a difference
M1.7 <- blmer(DHA_perc ~ Mode*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.7) ## See correlation of fixed effects
sjPlot::tab_model(M1.7)
BIC(M1.7)
car::Anova(M1.7,"III")
performance::check_collinearity(M1.7) ## VIF > 4

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Mode | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Country | Mode))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(DHA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(DHA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## See influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),],main="Belgium")
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]))
## Norway
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),],main="Norway")
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]))
## Sweden
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),],main="Sweden")
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]))
## Romania
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),],main="Romania")
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.7,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA:ALA - CORRELATION WITH KEY TRAITS                          ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_ALA_ratio ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_ALA_ratio ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]))

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
M1.5b <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std[-c(47,152,181),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)-BIC(M1.5b) ## Considerable improvement
plot(M1.5) ## See three extreme values
which.min(resid(M1.5))
which.max(resid(M1.5))

## Post-hoc test: Country*Site
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_ALA_ratio ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##          MODE          ##
############################

## Test if Mode makes a difference
M1.7 <- blmer(EPA_ALA_ratio ~ Mode*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.7) ## See correlation of fixed effects
sjPlot::tab_model(M1.7)
BIC(M1.7)
car::Anova(M1.7,"III")
performance::check_collinearity(M1.7) ## VIF > 4

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Mode | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.7, ~ Country | Mode))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_ALA_ratio ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(EPA_ALA_ratio ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## See influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_std[which(Spider_LME_data_std$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.7,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                               VARIATION PARTITIONING                                    ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##             VARIATION PARTITIONING - EXPLORE INFLUENCE OF COUNTRY                       ##        
##                                                                                         ##    
#############################################################################################

##*****************
head(Spider_LME_data_std)

Connectivity <- Spider_LME_data_std$EPA_perc

##*****************
## Connectivity
varpart_mod1 <- varpart(Connectivity,
                        ~ Spider_LME_data_std$dis4,
                        ~ Spider_LME_data_std$Country,
                        ~ Spider_LME_data_std$Genus_family)
showvarparts(3)
plot(varpart_mod1)
varpart_mod1

Shared_variation <- 0.11345 + -0.00060 + 0.13936 
Shared_variation 
## 0.25221 = 25.22%

Gmm1 <- model.matrix( ~ dis4, Spider_LME_data_std) 
Gmm2 <- model.matrix( ~ Country, Spider_LME_data_std) 
Gmm3 <- model.matrix( ~ Genus_family, Spider_LME_data_std) 

rda_result <- rda(Connectivity ~ Gmm1)
RsquareAdj(rda_result) ## 0.2505548

rda_result <- rda(Connectivity ~ Gmm1 + Gmm2 + Gmm3)
RsquareAdj(rda_result)
vif.cca(rda_result)
## 0.6169633

rda_result <- rda(Connectivity ~ Gmm1 + Condition(Gmm2) + Condition(Gmm3))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0 NS

rda_result <- rda(Connectivity ~ Gmm2 + Condition(Gmm1) + Condition(Gmm3))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0.1491029

rda_result <- rda(Connectivity ~ Gmm3 + Condition(Gmm1) + Condition(Gmm2))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0.1382069

##*****************
## Connectivity
varpart_mod1 <- varpart(Connectivity,
                        ~ Spider_LME_data_std$dis3,
                        ~ Spider_LME_data_std$Country,
                        ~ Spider_LME_data_std$Genus_family)
showvarparts(3)
plot(varpart_mod1)
varpart_mod1

Shared_variation <- 0.06705 + -0.00205 + 0.05863 
Shared_variation 
## 0.12363 = 12.36%

Gmm1 <- model.matrix( ~ dis3, Spider_LME_data_std) 
Gmm2 <- model.matrix( ~ Country, Spider_LME_data_std) 
Gmm3 <- model.matrix( ~ Genus_family, Spider_LME_data_std) 

rda_result <- rda(Connectivity ~ Gmm1)
RsquareAdj(rda_result) ## 0.1272642

rda_result <- rda(Connectivity ~ Gmm1 + Gmm2 + Gmm3)
RsquareAdj(rda_result)
vif.cca(rda_result)
## 0.6222593

rda_result <- rda(Connectivity ~ Gmm1 + Condition(Gmm2) + Condition(Gmm3))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0.003635108 NS

rda_result <- rda(Connectivity ~ Gmm2 + Condition(Gmm1) + Condition(Gmm3))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0.1955113

rda_result <- rda(Connectivity ~ Gmm3 + Condition(Gmm1) + Condition(Gmm2))
anova(rda_result, step=999, perm.max=999)
RsquareAdj(rda_result)
## adj.r.sqr 0.1396563


#############################################################################################
##                                                                                         ##
##                                  WEB vs GROUND SPIDERS                                  ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                                       WEB SPIDERS                                       ##        
##                                                                                         ##    
#############################################################################################

## Relabel the Forest sites 
Spider_LME_data_web_std[which(Spider_LME_data_web_std$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_web_std[which(Spider_LME_data_web_std$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_web_std[which(Spider_LME_data_web_std$Type_1=="LDS"),6]   <- "FOR"

head(Spider_LME_data_web_std)

#############################################################################################
##                                                                                         ##
##                          ALA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(ALA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(ALA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(ALA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),],main="Belgium")
abline(lm(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]))
## Norway
plot(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),],main="Norway")
abline(lm(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]))
## Sweden
plot(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),],main="Sweden")
abline(lm(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]))
## Romania
plot(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),],main="Romania")
abline(lm(ALA_perc ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(ALA_perc ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]))

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(ALA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis3 and Country > 2

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(ALA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## See influence of interation betwwn SITE and COUNTRY
M1.6 <- blmer(ALA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## See difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(ALA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Mode makes a difference
M1.10 <- blmer(ALA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_web_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## No interaction

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.7,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
plot(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.1)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country < 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)
plot(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## Almost a difference
performance::check_collinearity(M1.6) ## VIF is ok
plot(M1.6)

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
plot(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
plot(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Body Size makes a difference
M1.10 <- blmer(EPA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_web_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## Almost an influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]))
## Norway
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),],main="Norway")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]))
## Sweden
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]))
## Romania
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),],main="Romania")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          DHA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
M1.1b <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)-BIC(M1.1b)
which.max(resid(M1.1))
plot(M1.1) ## See outlier
plot(M1.1b) ## See outlier
sjPlot::tab_model(M1.1b)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(DHA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(DHA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(DHA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
plot(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(DHA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5) ## Some more extreme values
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(DHA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
plot(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_web_std[-c(43),])
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_web_std[-c(43),]),col="red")


## Test if Body Size interacts with Country
M1.9 <- blmer(DHA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std[-c(43),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(DHA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_web_std[-c(43),],
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## No influence of interaction

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA:ALA - CORRELATION WITH KEY TRAITS                          ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
plot(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_ALA_ratio ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]))

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_ALA_ratio ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ dis3, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis3 + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),]))

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
M1.5b <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
               data=Spider_LME_data_web_std[-c(25),],
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)-BIC(M1.5b) ## Considerable improvement
plot(M1.5b) ## See three extreme values
which.min(resid(M1.5))
which.max(resid(M1.5))

## Post-hoc test: Country*Site
L.S <- pairs(lsmeans(M1.5b, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5b, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_ALA_ratio ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_web_std[-c(25),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5b,M1.6) ## There is a difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Type_1))
test(L.S, adjust = "tukey")

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std[-c(25),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
plot(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_ALA_ratio ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_web_std[-c(25),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
plot(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Mode makes a difference
M1.10 <- blmer(EPA_ALA_ratio ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_web_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
plot(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## See influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_web_std[which(Spider_LME_data_web_std$Country=="RO"),],main="Romania")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                                       GROUND SPIDERS                                    ##        
##                                                                                         ##    
#############################################################################################

## Relabel the Forest sites 
Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Type_1=="LDS"),6]   <- "FOR"

head(Spider_LME_data_ground_std)

#############################################################################################
##                                                                                         ##
##                          ALA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(ALA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
plot(M1.1)
which.max(resid(M1.1))
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(ALA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(ALA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(ALA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis3 and Country > 2

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(ALA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## See influence of interation betwwn SITE and COUNTRY
M1.6 <- blmer(ALA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## See difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(ALA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
plot(M1.8)
car::Anova(M1.8,"II") ## Almost significant

## Test if Body Size interacts with Country
M1.9 <- blmer(ALA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-c(8,91,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Mode makes a difference
M1.10 <- blmer(ALA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_ground_std[-c(8,91,93),],
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## No interaction

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
plot(M1.1)
BIC(M1.1)
plot(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.1)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country < 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)
plot(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## Almost a difference
performance::check_collinearity(M1.6) ## VIF is ok
plot(M1.6)

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
plot(M1.8)
car::Anova(M1.8,"II")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
plot(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Body Size makes a difference
M1.10 <- blmer(EPA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_ground_std,
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF > 4
anova(M1.9,M1.10) ## Almost an influence of interaction

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="BE"),],main="Belgium")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="BE"),]))
## Norway
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="NO"),],main="Norway")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="NO"),]))
## Sweden
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="SE"),],main="Sweden")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="SE"),]))
## Romania
plot(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="RO"),],main="Romania")
abline(lm(EPA_perc ~ Body_Size, data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_perc ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std[which(Spider_LME_data_ground_std$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          DHA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
M1.1b <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)-BIC(M1.1b)
which.max(resid(M1.1))
plot(M1.1) ## See outlier
plot(M1.1b) ## See outlier
sjPlot::tab_model(M1.1b)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(DHA_perc ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(DHA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(DHA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
plot(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(DHA_perc ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5) ## Some more extreme values
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(DHA_perc ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(DHA_perc ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
plot(M1.8)
BIC(M1.8)
car::Anova(M1.8,"II")
plot(DHA_perc ~ Body_Size, data=Spider_LME_data_ground_std[-c(31,74),])
abline(lm(DHA_perc ~ Body_Size, data=Spider_LME_data_ground_std[-c(31,74),]),col="red")

## Test if Body Size interacts with Country
M1.9 <- blmer(DHA_perc ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-c(31,74),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 4

## Test if Mode makes a difference
M1.10 <- blmer(DHA_perc ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_ground_std[-c(31,74),],
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
summary(M1.10)
BIC(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## No influence of interaction

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA:ALA - CORRELATION WITH KEY TRAITS                          ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(81,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
which.max(resid(M1.1))
which.min(resid(M1.1))
plot(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2
Spider_LME_data_ground_std[c(81,93),c(1:10)]

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(81,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_ALA_ratio ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(81,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_ALA_ratio ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std[-c(81,93),],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
plot(M1.5) ## See three extreme values
which.min(resid(M1.5))
which.max(resid(M1.5))
performance::check_collinearity(M1.5) ## Check VIF of without interaction

## Post-hoc test: Country*Site
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_ALA_ratio ~ Type_1 + Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_ground_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## There is a difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Type_1))
test(L.S, adjust = "tukey")

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

############################
##       BODY SIZE        ##
############################

## Test influence of Body Size
M1.8 <- blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-81,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.8) ## See correlation of fixed effects
sjPlot::tab_model(M1.8)
BIC(M1.8)
plot(M1.8)
which.max(resid(M1.8)) ## -81
car::Anova(M1.8,"II")

## Trophic connectivity increases with spider body size
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std[-81,])
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std[-81,]),col="red")

## Test if Body Size interacts with Country
M1.9 <- blmer(EPA_ALA_ratio ~ Body_Size*Country + (1|Site_block/Site_name), 
              data=Spider_LME_data_ground_std[-81,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
summary(M1.9) ## See correlation of fixed effects
sjPlot::tab_model(M1.9)
BIC(M1.9)
plot(M1.9)
car::Anova(M1.9,"III")
performance::check_collinearity(M1.9) ## VIF > 3

## Test if Body Size and Country makes a difference
M1.10 <- blmer(EPA_ALA_ratio ~ Body_Size + Country + (1|Site_block/Site_name), 
               data=Spider_LME_data_ground_std[-81,],
               lmerControl(optimizer = "Nelder_Mead"),
               REML=T)
sjPlot::tab_model(M1.10)
BIC(M1.10)
plot(M1.10)
car::Anova(M1.10,"II")
performance::check_collinearity(M1.10) ## VIF < 2
anova(M1.9,M1.10) ## See influence of interaction

## Create data set without outlier
Spider_LME_data_ground_std_mOUT <- Spider_LME_data_ground_std[-81,]

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ Body_Size, data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ Body_Size + (1|Site_block/Site_name), data=Spider_LME_data_ground_std_mOUT[which(Spider_LME_data_ground_std_mOUT$Country=="RO"),]))

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                                      LYCOSIDAE SPIDERS                                  ##        
##                                                                                         ##    
#############################################################################################

## Relabel the Forest sites 
Spider_LME_data_Lycos_std[which(Spider_LME_data_Lycos_std$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_Lycos_std[which(Spider_LME_data_Lycos_std$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_Lycos_std[which(Spider_LME_data_Lycos_std$Type_1=="LDS"),6]   <- "FOR"

head(Spider_LME_data_Lycos_std)

#############################################################################################
##                                                                                         ##
##                          ALA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(ALA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
plot(M1.1)
which.max(resid(M1.1))
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(ALA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 2 - raises a red flag
car::Anova(M1.2,"II")

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(ALA_perc ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(ALA_perc ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis3 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(ALA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## See influence of interation betwwn SITE and COUNTRY
M1.6 <- blmer(ALA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## See difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
plot(M1.1)
BIC(M1.1)
plot(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.1)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_perc ~ dis4*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_perc ~ dis3*Country + (1|Site_block/Site_name) + (1|Genus_family), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country < 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)
plot(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## There is no interaction
performance::check_collinearity(M1.6) ## VIF is ok
plot(M1.6)

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

#############################################################################################
##                                                                                         ##
##                          DHA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
which.max(resid(M1.1))
plot(M1.1) ## See outlier
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(DHA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(DHA_perc ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(DHA_perc ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
plot(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(DHA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5) ## Some more extreme values
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(DHA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std[-44,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

#############################################################################################
##                                                                                         ##
##                          EPA:ALA - CORRELATION WITH KEY TRAITS                          ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
plot(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_ALA_ratio ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_ALA_ratio ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
plot(M1.5) ## See three extreme values
performance::check_collinearity(M1.5) ## Check VIF without interaction

## Post-hoc test: Country*Site
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_ALA_ratio ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Lycos_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## There is a difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Type_1))
test(L.S, adjust = "tukey")

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

#############################################################################################
##                                                                                         ##
##                                 TETRAGNATHIDAE SPIDERS                                  ##        
##                                                                                         ##    
#############################################################################################

## Relabel the Forest sites 
Spider_LME_data_Tetra_std[which(Spider_LME_data_Tetra_std$Type_1=="FOR1"),6]  <- "FOR"
Spider_LME_data_Tetra_std[which(Spider_LME_data_Tetra_std$Type_1=="FOR2"),6]  <- "FOR"
Spider_LME_data_Tetra_std[which(Spider_LME_data_Tetra_std$Type_1=="LDS"),6]   <- "FOR"

head(Spider_LME_data_Tetra_std)

#############################################################################################
##                                                                                         ##
##                          ALA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(ALA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
plot(M1.1)
which.max(resid(M1.1))
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(ALA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for Country just > 2 - raises a red flag
car::Anova(M1.2,"II")

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(ALA_perc ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(ALA_perc ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis3 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(ALA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## See influence of interation betwwn SITE and COUNTRY
M1.6 <- blmer(ALA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## See difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6,M1.8,M1.9,M1.10)

#############################################################################################
##                                                                                         ##
##                          EPA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
plot(M1.1)
BIC(M1.1)
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.1)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country (just) > 2 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_perc ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_perc ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country < 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
BIC(M1.5)
plot(M1.5)

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## Marginally no interaction
performance::check_collinearity(M1.6) ## VIF is ok
plot(M1.6)

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

#############################################################################################
##                                                                                         ##
##                          DHA - CORRELATION WITH KEY TRAITS                              ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(DHA_perc ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
BIC(M1.1)
which.max(resid(M1.1))
plot(M1.1) ## See outlier
performance::check_collinearity(M1.1) ## VIF for dis4 and dis3 < 2
Spider_LME_data_Tetra_std[19,c(1:3,5:8,33)] ## Extreme DHA value SKA
head(Spider_LME_data_Tetra_std)

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(DHA_perc ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
plot(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(DHA_perc ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(DHA_perc ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
plot(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(DHA_perc ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
plot(M1.5) ## Some more extreme values
car::Anova(M1.5,"III")
BIC(M1.5)

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(DHA_perc ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-19,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No difference
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

#############################################################################################
##                                                                                         ##
##                          EPA:ALA - CORRELATION WITH KEY TRAITS                          ##        
##                                                                                         ##    
#############################################################################################

## MAIN MODEL - see methods for justification of steps to deal with colinearity issues
## Test the correlation of dis3 and dis4
M1.1 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Tetra_std,
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
M1.1b <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.1)
sjPlot::tab_model(M1.1b) ## Considerable improvement, suggests genuine outlier
plot(M1.1)
plot(M1.1b)
BIC(M1.1)-BIC(M1.1b)
which.min(resid(M1.1))
performance::check_collinearity(M1.1b) ## VIF for dis4 and dis3 < 2

## Test the correlation of dis3 and dis4 with Country added as a covariate
M1.2 <- blmer(EPA_ALA_ratio ~ dis3 + dis4 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.2)
BIC(M1.2)
summary(M1.2)
performance::check_collinearity(M1.2) ## VIF for dis4 and Country > 3 - raises a red flag

## Test the correlation of dis4 with Country added as a interaction term
M1.3 <- blmer(EPA_ALA_ratio ~ dis4*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.3)
BIC(M1.3)
performance::check_collinearity(M1.3) ## VIF for dis4:Country and Country > 10
car::Anova(M1.3,"III")

Spider_LME_data_Tetra_std_mOUT <- Spider_LME_data_Tetra_std[-3,]

## Plot for interpretation
par(mfrow=c(2,2))
## Belgium
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="BE"),],main="Belgium")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="BE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block), data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="BE"),]))
## Norway
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="NO"),],main="Norway")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="NO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block), data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="NO"),]))
## Sweden
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="SE"),],main="Sweden")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="SE"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block), data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="SE"),]))
## Romania
plot(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="RO"),],main="Romania")
abline(lm(EPA_ALA_ratio ~ dis4, data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="RO"),]),col="red")
sjPlot::tab_model(blmer(EPA_ALA_ratio ~ dis4 + (1|Site_block), data=Spider_LME_data_Tetra_std_mOUT[which(Spider_LME_data_Tetra_std_mOUT$Country=="RO"),]))

## Assess individual correlations for each trait
## Test the correlation of dis3 with Country added as a interaction term
M1.4 <- blmer(EPA_ALA_ratio ~ dis3*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.4)
BIC(M1.4)
performance::check_collinearity(M1.4) ## dis4 and Country > 2
car::Anova(M1.4,"III")

############################
##        SITE TYPE       ##
############################

## Assess individual correlations for each trait
## First test the difference betwween site types
M1.5 <- blmer(EPA_ALA_ratio ~ Type_1*Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.5)
car::Anova(M1.5,"III")
plot(M1.5) ## See three extreme values
performance::check_collinearity(M1.5) ## Check VIF without interaction

## Post-hoc test: Country*Site
L.S <- pairs(lsmeans(M1.5, ~ Type_1 | Country))
test(L.S, adjust = "tukey")

## Post-hoc test: Country*Mode
L.S <- pairs(lsmeans(M1.5, ~ Country | Type_1))
test(L.S, adjust = "tukey")

## No significant interaction, but see influence of COUNTRY
M1.6 <- blmer(EPA_ALA_ratio ~ Type_1 + Country + (1|Site_block), 
              data=Spider_LME_data_Tetra_std[-3,],
              lmerControl(optimizer = "Nelder_Mead"),
              REML=T)
sjPlot::tab_model(M1.6)
car::Anova(M1.6,"II")
BIC(M1.6)
anova(M1.5,M1.6) ## No interaction 
performance::check_collinearity(M1.6) ## VIF is ok

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Type_1))
test(L.S, adjust = "tukey")

## Post-hoc test: Country
L.S <- pairs(lsmeans(M1.6, ~ Country))
test(L.S, adjust = "tukey")

## Tidyup workspace
rm(M1.1,M1.2,M1.3,M1.4,M1.5,M1.6)

