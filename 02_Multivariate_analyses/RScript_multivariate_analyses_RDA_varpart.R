###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Figs.4 and 5 in "Fatty acid biomarkers reveal landscape influences on         ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Make plots for ground and web spiders, and Lycosidae and Tetragnathidae
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 4 July 2021
# Updated: 7 June 2025

## This R script performs analyses and creates figures relating to multivariate analyses of spider fatty acid composition
## Includes redundancy analyses (RDA) and variation partitioning (varpart) - See Figures 4 and 5
## Also includes permutational multivariate analysis of variance (PERMANOVA)

## Metadata for "Burdon_et_al_Figs.4_5_multivariate_data.csv" (Numbers = column)
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
## 13-49. FA1-FA37:	Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
## 50.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 51.    chl_a_m2_day:	Algal accrual on tiles
## 52.    Width:	Stream wetted channel width (m)
## 53.    Riparian_PC1:	Riparian condition PC1 (36.7%)
## 54.    Catchment_PC1: Catchment impact PC1 (42.3%)	
## 55-87. s1-wnb5:	CWM abundances for 33 invertebrate dispersal trait modalities (e.g., aerial active, body size)
## 88.    Body_Size:	Mean body size of spider Family sampled for FAs

# Metadata for "Burdon_et_al_Figs.4_5_FA_lookup_table.csv" (Numbers = column)
## 1.     FA: Code for fatty acids (FA1-FA38)
## 2.     Short_name_2: Common chemical nomenclature used for fatty acids (FAs)
## 3.     FA_1: FA grouping 1 (separates out omega-3 and omega-6 PUFAs from other FA groups)
## 4.     FA_2: FA grouping 2 (separates out PUFAs from other FA groups)
## 5.     Resource: Resources typically associated with FA
## 6.     Specific_PUFA: Short name for selected PUFAs
## 7.     Long_name: long name of FA

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
library(rJava)
library(venneuler)
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
library(packfor)

#############################################################################################
##                                                                                         ##
##                                       LOAD DATA                                         ##        
##                                                                                         ##    
#############################################################################################

## Load data for analyses
Spider_FA_data_new <- read_csv("Burdon_et_al_Figs.4_5_multivariate_data.csv")

## Load FA lookup for names of fatty acids 
FA_lookup <- read_csv("Burdon_et_al_Figs.4_5_FA_lookup_table.csv")

#############################################################################################
##                                                                                         ##
##                      REDUNDANCY ANALYSIS AND VARIATION PARTITIONING                     ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                            TRANSFORMATIONS AND STANDARDISATION                          ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset
Spider_FA_data_log <- Spider_FA_data_new

## Transform FAs by taking square-root of relative abundances (Hellinger)
Spider_FA_data_log[,c(13:49)] <- sqrt(Spider_FA_data_log[,c(13:49)])

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- na.omit(Spider_FA_data_log[,-12])

## Remove NAs (256 -> 238 samples)
Spider_FA_data_log <- na.omit(Spider_FA_data_log)

## Check the number of sites used
Data_sites <- Spider_FA_data_log %>%
  group_by(Country)%>%
  count(Site_name)

## Create data frame and add one for each site
Data_sites$Site <- 1

## Count sites from each Country
Data_sites %>%
  group_by(Country)%>%
  count(Site)

## Transformations predictor variables where appropriate
Spider_FA_data_log[,c(50)] <- log1p(Spider_FA_data_log[,c(50)]) ## Log+1 transform chl-a concentrations per day
Spider_FA_data_log[,c(51)] <- log(Spider_FA_data_log[,c(51)]) ## Log transform stream wetted width (m)
Spider_FA_data_log[,c(87)] <- log(Spider_FA_data_log[,c(87)]) ## Log-transform body size

## Standardize predictors
Spider_FA_data_log[,c(49:87)] <- decostand(Spider_FA_data_log[,c(49:87)],"standardize",na.action=na.exclude)

## Check output
colnames(Spider_FA_data_log)

## Check correlations using hypothesized drivers
head(Spider_FA_data_log[,c(49:87)])
global <- Spider_FA_data_log[,c(49:53,56,61,64,68,69,70,71,83,84,86,87)]
head(global)
myMat.cor <- cor(global, method=c("pearson"))

## Note highly correlated traits:
## s3 ~ cd1 > 0.8
## s3 ~ life1 > 0.8
## life1 ~ cd1 > 0.8
## wnb2 ~ dis3 > 0.8

## Need to drop wnb2 ##83 (Diptera) and 
## s3 ##56 (Ephemeroptera, Trichoptera) and 
## life1 ##70
## cy2 and cd1
## NOTE: Also drop dis4 since VIF > 4
Gmm <- model.matrix( ~ Spider_body_size + chl_a_m2_day + Width + Catchment_PC1 + Riparian_PC1 + life2 + dis3 + wnb3 + wnb5 + Body_Size, Spider_FA_data_log) 
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.2513519

## Subset data
head(Spider_FA_data_log[,c(49:87)])
global <- Spider_FA_data_log[,c(49:53,68,71,84,86,87)]
head(global)

#calculate the pearson correlation coefficient matrix
## Drop metals due to NAs
myMat.cor <- cor(global, method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)

## Blanchet et al. (2008)
## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Spider_FA_data_log[,c(12:48)], global, nperm=999, adjR2thresh = 0.2513519, alpha = 0.05)
as.matrix(M1$variables)
#rm(M1,global)
## "Catchment_PC1"   
## "Body_Size"       
## "Riparian_PC1"    
## "wnb5"            
## "chl_a_m2_day"    
## "wnb3"            
## "dis3"            
## "Spider_body_size"

## Note: p-value for "life2" frequently between 0.052 and 0.054, so should consider if it adds value

## OLD - UPDATED DATA
## "cy2"             
## "Riparian_PC1"    
## "Body_Size"       
## "wnb5"            
## "wnb3"            
## "cd1"             
## "chl_a_m2_day"    
## "Catchment_PC1"   
## "life2"           
## "Spider_body_size"

## Check model
Gmm <- model.matrix( ~ Spider_body_size + Body_Size + chl_a_m2_day + Riparian_PC1 + Catchment_PC1 + dis3 + wnb3 + wnb5, Spider_FA_data_log) 
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
extractAIC(rda_result)
## adj.r.sqr 0.2461254

## head(Spider_FA_data_log)
plot_data <- Spider_FA_data_log[,c(49,50,52,53,68,84,71,86,87)]
colnames(plot_data)
names(plot_data) <- c("SpiderCWM","Algae","Riparian","Catchment","dis3","wnb3","life2","wnb5","SpiderBS") 

## Create data frame for PUFA data
dSpider <- Spider_FA_data_log[,c(12:48)]

# Analyse and plot the RDA using ggplot (ggord package)
## Drop cy2 VIF > 4 (Not relevant in revised model since input data has been updated)
rda_result <- rda(dSpider ~ SpiderCWM + Algae + Riparian + Catchment + life2 + dis3 + wnb3 + wnb5 + SpiderBS, plot_data)
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
anova.cca(rda_result, by="axis")
vif.cca(rda_result)
extractAIC(rda_result)
##    10.000 -872.4041 (with life2)
##      9.00 -872.4015 (without life2)
## Change in AIC -872.4015--872.4041 = 0.0026 (some support for keep life2 in)

## Create objects from RDA
RDA_rsq_result <- RsquareAdj(rda_result)
RDA_permtest_result <- anova(rda_result, step=999, perm.max=999)
RDA_axes_result <- anova.cca(rda_result, by="axis")
RDA_vif_result <- vif.cca(rda_result)

#############################################################################################
##                                                                                         ##
##                               VISUALISE THE RDA RESULTS                                 ##                
##                                                                                         ##    
#############################################################################################

## Test data
sites <- scores(rda_result , display = "sites", choices=c(1,2))
species <- scores(rda_result , display = "species", choices=c(1,2))
vectors <- scores(rda_result , display = "bp", choices=c(1,2))
RDA_sites <- data.frame(Spider_FA_data_log[,c(1:8)], sites)

RDA_FA <- as.data.frame(species)## Create dataframe
RDA_FA$FA <- rownames(species)## Create new taxa column
head(RDA_FA)
RDA_FA <- merge(RDA_FA,FA_lookup[,c(1:2)],by="FA")
RDA_FA <- RDA_FA[,c(4,2,3)]

## Take the column sums of relative abundance data
rowSums(Spider_FA_data_new[,c(13:49)])

## Calculate the mean abundance of the FAs
mean_abundance <- as.data.frame(colSums(Spider_FA_data_new[,c(13:49)]*100)/256)

## Change row names to a column
mean_abundance$FA <- row.names(mean_abundance) 

## Rename the first column
colnames(mean_abundance)[1] <- "Mean_percentage"

## Sort highest to lowest
sorted <- mean_abundance[order(-mean_abundance$Mean_percentage),]

## Sum for 95%
sum(sorted$Mean_percentage[1:10])
FA_sorted <- merge(sorted,as.data.frame(FA_lookup[,c(1,2)]),by="FA",sort = FALSE)
#FA_sorted <- FA_sorted[order(-FA_sorted$Mean_percentage),]
unique(FA_sorted$Short_name_2[1:10])

Influential_FA  <-  c("C18.1n-9","C16.0","C18.0","C18.2n-6","C14.0","C20.5n-3","C16.1n-7","C17.0","C20.4n-6","C20.1n-9")

## Top 10 FAs in terms of mean abundance 
FA_influential <- RDA_FA[RDA_FA$Short_name_2 %in% Influential_FA,]

sites <- scores(rda_result, display = "sites", choices=c(1,2))
scrs <- data.frame(as.data.frame(sites), Site_block = Spider_FA_data_log$Site_block)
cent <- aggregate(cbind(RDA1, RDA2) ~ Site_block, data = scrs, FUN = mean)
scrs <- scrs[,c(3,1,2)]
scrs <- setNames(scrs, c("Site_block","RDA1","RDA2"))
cent <- setNames(cent, c("Site_block","oRDA1","oRDA2"))
segs <- merge(scrs, cent, by.x = "Site_block", sort = FALSE)

##***********************************************  
##      CREATE PLOTS TO VISUALIZE RESULTS       
##***********************************************

##***********************************************  
##       RDA: PLOT SPIDER %FA COMPOSITION      
##***********************************************

## Create plot for spider FA composition with 10 most abundant FAs shown
RDA_plot <- ggplot(data=scrs, aes(RDA1, RDA2)) + 
  # Plot equivalent of ordispider
  #geom_segment(data = segs, mapping = aes(xend = oRDA1, yend = oRDA2),linetype=1, color="grey50") +
  # Plot points
  geom_point(data=RDA_sites, aes(RDA1, RDA2, colour=Country, shape=Type_1), size=4, alpha=0.8) + 
  # Plot ellipses with fill
  stat_ellipse(data=RDA_sites, aes(RDA1, RDA2, colour=Country, fill = Country),geom = "polygon", alpha = 0.25) +
  # Create shape for site types
  scale_shape_manual(breaks=c("UBF", "FBF", "FOR"),
                     values=c(15, 17, 16),
                     labels = c("Unbuffered", "Buffered", "Forest"),
                     name = "Site type"
                     #, guide = FALSE
  ) +
  # Create colours for points
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"),
                      name = "Country"
                      #, guide = FALSE
  ) +
  # Create fill colours for the ellipses
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    name = "Country"
                    #, guide = FALSE
  ) +
  # Create axes labels
  labs(x="RDA1 (52.1%)", y="RDA2 (22.2%)") +
  # Expand x-axis label
  coord_cartesian(xlim = c(-0.32, 0.35)) +
  # Plot the 10 most abundant FAs
  geom_text(data=FA_influential[-7,], aes(x=RDA1, y=RDA2, label=Short_name_2), color="blue", size=4.5) +
  geom_text(data=FA_influential[7,], aes(x=RDA1-0.02, y=RDA2, label=Short_name_2), color="blue", fontface="bold", size=4.5) +
  # Clear out backgrounds and gridlines
  theme_bw() +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=13),
        axis.title.y=element_text(vjust=1.25, size = 18, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 18, face="bold", color="black"),
        axis.text.y=element_text(size = 14, color="black"),
        axis.text.x=element_text(size = 14, color="black"))

##***********************************************  
##             RDA: PLOT PREDICTORS     
##***********************************************

RDA_vectors <- as.data.frame(vectors)## Create dataframe
RDA_vectors$vectors <- rownames(RDA_vectors)## Create new taxa column

p_Vectors <- ggplot() + 
  
  geom_hline(yintercept=0, linetype="dotted", color = "grey20") +
  
  geom_vline(xintercept=0, linetype="dotted", color = "grey20") +
  
  geom_segment(data=vectors, aes(x=0, y=0, xend=RDA1, yend=RDA2), arrow=arrow(length=unit(0.3,"cm"))) +
  
  geom_text(data=RDA_vectors[1,], aes(x=RDA1-0.1, y=RDA2-0.05, label=vectors), size=4.5) + ## SpiderCWM
  geom_text(data=RDA_vectors[2,], aes(x=RDA1, y=RDA2+0.1, label=vectors), size=4.5) + ## Algae
  geom_text(data=RDA_vectors[3,], aes(x=RDA1-0.25, y=RDA2, label=vectors), size=4.5) + ## Riparian
  geom_text(data=RDA_vectors[4,], aes(x=RDA1+0.18, y=RDA2+0.18, label=vectors), size=4.5) + ## Catchment
  geom_text(data=RDA_vectors[5,], aes(x=RDA1-0.05, y=RDA2+0.1, label=vectors), size=4.5) + ## life2
  geom_text(data=RDA_vectors[9,], aes(x=RDA1+0.15, y=RDA2-0.05, label=vectors), size=4.5) + ## SpiderBS
  geom_text(data=RDA_vectors[6,], aes(x=RDA1, y=RDA2-0.075, label=vectors), size=4.5) + ## dis3
  geom_text(data=RDA_vectors[7,], aes(x=RDA1, y=RDA2-0.075, label=vectors), size=4.5) + ## wnb3
  geom_text(data=RDA_vectors[8,], aes(x=RDA1-0.05, y=RDA2-0.075, label=vectors), size=4.5) + ## wnb5
  
  #geom_text(data=FA_influential[7,], aes(x=RDA1, y=RDA2, label=Short_name_2), color="blue") +
  
  ## 
  labs(x="RDA1 (52.1%)", y="RDA2 (22.2%)") +
  
  theme_bw() +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 18, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 18, face="bold", color="black"),
        axis.text.y=element_text(size = 14, color="black"),
        axis.text.x=element_text(size = 14, color="black"))

##***********************************************  
##          RDA: CREATE FINAL PLOT    
##***********************************************

## Create multiple
png(filename="Burdon_et_al_Fig.4.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=5, 
    pointsize=20, 
    res=600)

ggarrange(ggMarginal(RDA_plot, groupColour = TRUE, groupFill = TRUE), p_Vectors, widths = c(2, 1.4),
          ncol = 2, nrow = 1,
          labels = c("(a)", "(b)"), font.label = list(size = 18), label.x = -0.018)

dev.off()

#############################################################################################
##                                                                                         ##
##                                PLOT OTHER DIMENSIONS                                    ##                
##                                                                                         ##    
#############################################################################################

## Create plot data
plot_data <- Spider_FA_data_log[,c(49,50,52,53,68,84,71,86,87)]
colnames(plot_data)
names(plot_data) <- c("SpiderCWM","Algae","Riparian","Catchment","dis3","wnb3","life2","wnb5","SpiderBS") 

## Analyse and plot the RDA using ggplot (ggord package)
dSpider <- Spider_FA_data_log[,c(12:48)]
colnames(dSpider) <- c(FA_lookup$Short_name_2[-38])
rda_result <- rda(dSpider ~ SpiderCWM + Algae + Riparian + Catchment + life2 + dis3 + wnb3 + wnb5 + SpiderBS, plot_data)

## Plor figure
png(filename="Burdon_et_al_Fig.S2.png", 
    type="cairo",
    units="in", 
    width=9, 
    height=8, 
    pointsize=20, 
    res=600)

p1 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "3"), ellipse = TRUE,
            ellipse_pro = 0.95, ptslab=T, obslab=F, size= 5, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium")) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), guide = FALSE) +
  guides(fill=guide_legend(title="Country")) + 
  expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

p2 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("2", "3"), ellipse = TRUE,
            ellipse_pro = 0.95, ptslab=T, obslab=F, size= 5, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium")) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), guide = FALSE) +
  guides(fill=guide_legend(title="Country")) + 
  expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))


p3 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "4"), ellipse = TRUE,
            ellipse_pro = 0.95, ptslab=T, obslab=F, size= 5, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium")) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), guide = FALSE) +
  guides(fill=guide_legend(title="Country")) + 
  expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

p4 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("2", "4"), ellipse = TRUE,
            ellipse_pro = 0.95, ptslab=T, obslab=F, size= 5, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium")) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), guide = FALSE) +
  guides(fill=guide_legend(title="Country")) + 
  expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

p5 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("3", "4"), ellipse = TRUE,
            ellipse_pro = 0.95, ptslab=T, obslab=F, size= 5, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium")) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), guide = FALSE) +
  guides(fill=guide_legend(title="Country")) + 
  expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

## Arrange using grid plot
library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3, common.legend = TRUE, legend="right")

require(grid)
grid.text("(a)", x = unit(0.02, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(c)", x = unit(0.02, "npc"), y = unit(0.65, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(e)", x = unit(0.02, "npc"), y = unit(0.31, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(b)", x = unit(0.47, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(d)", x = unit(0.47, "npc"), y = unit(0.65, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))

dev.off()

############################################################################################
##                                                                                        ##
##                               VARIATION PARTITIONING                                   ##
##                                                                                        ## 
############################################################################################
############################################################################################
##                                                                                        ##
##                                TEST FOUR HYPOTHESES                                    ##
##                                                                                        ## 
############################################################################################

## Tests four hypotheses: algae productivity, insect dispersal, spider community, and environment
varpart_mod1 <- varpart(dSpider,
                        ~ plot_data$SpiderCWM + plot_data$SpiderBS,
                        ~ plot_data$dis3 + plot_data$life2 + plot_data$wnb3 + plot_data$wnb5,
                        ~ plot_data$Algae,
                        ~ plot_data$Riparian + plot_data$Catchment)
showvarparts(4)
varpart_mod1

## Create final plot for combined approach
png(filename="Burdon_et_al_Fig.5a.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=16, 
    res=800)

MyVenn <- venneuler(c(Spider=0.04261, Insect=0.04712, Algae=0.02018, Environment=0.04218,
                      "Spider&Insect"=0.00005,
                      "Insect&Algae"=0.00327,
                      "Spider&Algae"=0,
                      "Spider&Environment"=0.00793,
                      "Environment&Insect"=0.04304,
                      "Environment&Algae"=0.00348,
                      "Spider&Insect&Environment"=0.00565,
                      "Spider&Insect&Algae"=0.00105,
                      "Environment&Insect&Algae"=0.03272,
                      "Spider&Algae&Environment"=0.00115,
                      "Spider&Insect&Algae&Environment"=0
))

MyVenn$labels <- c("","","","")

plot(MyVenn, col=c('#B15928','orange','#1B9E77',"#1F78B4"))

#plot(MyVenn, col=c('antiquewhite4','orange','green4','yellow4'))

dev.off()

## Now test the significance of the independent fractions
mm1 <- model.matrix( ~ SpiderCWM + SpiderBS, plot_data) 
mm2 <- model.matrix( ~ dis3 + life2 + wnb3 + wnb5, plot_data)
mm3 <- model.matrix( ~ Algae, plot_data) 
mm4 <- model.matrix( ~ Riparian + Catchment, plot_data) 

## Overall
rda_result <- rda(dSpider ~ mm1 + mm2 + mm3 + mm4)
M1_overall_rsq <- RsquareAdj(rda_result)
M1_overall_perm <- anova(rda_result, step=999, perm.max=999)

## Spider
rda_result <- rda(dSpider ~ mm1 +  Condition(mm2) + Condition(mm3) + Condition(mm4))
M1_spider_rsq <- RsquareAdj(rda_result)
M1_spider_perm <- anova(rda_result, step=999, perm.max=999)

## Insect
rda_result <- rda(dSpider ~ mm2 +  Condition(mm1) + Condition(mm3) + Condition(mm4))
M1_insect_rsq <- RsquareAdj(rda_result)
M1_insect_perm <- anova(rda_result, step=999, perm.max=999)

## Algae
rda_result <- rda(dSpider ~ mm3 +  Condition(mm2) + Condition(mm1) + Condition(mm4))
M1_algae_rsq <- RsquareAdj(rda_result)
M1_algae_perm <- anova(rda_result, step=999, perm.max=999)

## Environment
rda_result <- rda(dSpider ~ mm4 +  Condition(mm2) + Condition(mm3) + Condition(mm1))
M1_envr_rsq <- RsquareAdj(rda_result)
M1_envr_perm <- anova(rda_result, step=999, perm.max=999)

############################################################################################
##                                                                                        ##
##                        TEST INFLUENCE OF RIPARIAN CONDITION                            ##
##                                                                                        ## 
############################################################################################

## Variation partitioning - to test the independent contribution of riparian condition 
## Reduced approach
## Riparian, Food web (algae, insects, spider), Environment
varpart_mod2 <- varpart(dSpider,
                        ~ plot_data$Riparian ,
                        ~ plot_data$Algae + plot_data$dis3 + plot_data$life2 + plot_data$wnb3 + plot_data$wnb5 + plot_data$SpiderCWM + plot_data$SpiderBS,
                        ~ plot_data$Catchment)
showvarparts(3)
plot(varpart_mod2)
varpart_mod2

## NEWER VERSION: DIFFERENCE COLOUR SCHEME

## Create final plot for combined approach
png(filename="Burdon_et_al_Fig.5b.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=16, 
    res=800)

MyVenn <- venneuler(c(Riparian=0.01694, Foodweb=0.11407, Environment=0.01031,
                      "Riparian&Foodweb"=0.01470,
                      "Foodweb&Environment"=0.05899,
                      "Riparian&Environment"=0.01494,
                      "Riparian&Foodweb&Environment"=0.01923))

MyVenn$labels <- c("","","")

plot(MyVenn, col=c('#666666','peru','#7570B3'))

dev.off()

## Test data
mm1 <- model.matrix( ~ Riparian, plot_data) 
mm2 <- model.matrix( ~ Algae + dis3 + life2 + wnb3 + wnb5 + SpiderBS + SpiderCWM, plot_data)
mm3 <- model.matrix( ~ Catchment, plot_data) 

## Overall
rda_result <- rda(dSpider ~ mm1 + mm2 + mm3)
M1_overall_rsq <- RsquareAdj(rda_result)
M1_overall_perm <- anova(rda_result, step=999, perm.max=999)

## Riparian
rda_result <- rda(dSpider ~ mm1 +  Condition(mm2) + Condition(mm3))
M1_riparian_rsq <- RsquareAdj(rda_result)
M1_riparian_perm <- anova(rda_result, step=999, perm.max=999)

## Food web
rda_result <- rda(dSpider ~ mm2 +  Condition(mm1) + Condition(mm3))
M1_foodweb_rsq <- RsquareAdj(rda_result)
M1_foodweb_perm <- anova(rda_result, step=999, perm.max=999)

## Catchment
rda_result <- rda(dSpider ~ mm3 +  Condition(mm2) + Condition(mm1))
M1_envr_rsq <- RsquareAdj(rda_result)
M1_envr_perm <- anova(rda_result, step=999, perm.max=999)

############################################################################################
##                                                                                        ##
##                             Account for Spider Identity                                ##
##                                                                                        ## 
############################################################################################

varpart_mod2.2 <- varpart(dSpider,
                          ~ plot_data$Catchment + plot_data$Riparian + plot_data$Algae + 
                            plot_data$dis3 + plot_data$life2 + plot_data$wnb3 + plot_data$wnb5 + 
                            plot_data$SpiderCWM,
                          ~ Spider_FA_data_log$Genus_family)

plot(varpart_mod2.2)


## Test data
Gmm1 <- model.matrix( ~ Spider_body_size + dis3 + life2 + wnb3 + wnb5 + chl_a_m2_day + 
                        Riparian_PC1 + Catchment_PC1, Spider_FA_data_log) 
Gmm2 <- model.matrix( ~ Genus_family, Spider_FA_data_log) 
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 + Gmm2)
RsquareAdj(rda_result)
vif.cca(rda_result) ## Note Spiders that have high VIF values
## $adj.r.squared
## 0.322384

## Predictors
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 +  Condition(Gmm2))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.1549222

## Spider identity
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm2 +  Condition(Gmm1))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.1093458

## 0.322384 - (0.1549222 + 0.1093458)
## 0.058116

## 32.2% Combined influence of variation in PUFA 
## 15.5% Predictors independent of spider identity
## 10.9% Spider identity independent of predictors 
## 5.8% Shared influence


############################################################################################
##                                                                                        ##
##                        Account for Spider Identity and Country                         ##
##                                                                                        ## 
############################################################################################

## Need to drop Body Size (since differences in size match spider identifies)
varpart_mod3 <- varpart(dSpider,
                        ~ plot_data$Catchment + plot_data$Riparian + plot_data$Algae + 
                          plot_data$dis3 + plot_data$life2 + plot_data$wnb3 + plot_data$wnb5 + 
                          plot_data$SpiderCWM,
                        ~ Spider_FA_data_log$Genus_family,
                        ~ Spider_FA_data_log$Country)

plot(varpart_mod3)

## Include Country and Spider identity
## Account for Spider Identity - selected parameters
Gmm1 <- model.matrix( ~ Spider_body_size + dis3 + life2 + wnb3 + wnb5 + chl_a_m2_day + 
                        Riparian_PC1 + Catchment_PC1, Spider_FA_data_log) 
Gmm2 <- model.matrix( ~ Genus_family, Spider_FA_data_log)
Gmm3 <- model.matrix( ~ Country, Spider_FA_data_log) 
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 + Gmm2 + Gmm3)
RsquareAdj(rda_result)
vif.cca(rda_result)
## $adj.r.squared
## 0.4352103 ## Note High VIF values > 4 for Countries

## Predictors
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 +  Condition(Gmm2) + Condition(Gmm3))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.03027383

## Spider identity
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~  Gmm2 +  Condition(Gmm1) + Condition(Gmm3))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.1004845

## Country
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~  Gmm3 +  Condition(Gmm1) + Condition(Gmm2))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.1128263

## 0.4352103 - (0.03027383 + 0.1004845 + 0.1128263)
## 0.1916257

## 43.5% Combined influence of variation in PUFA 
## 3.0% Predictors independent of spider identity and country
## 10.0% Spider identity independent of environment and country
## 11.3% Country independent of environment and spider identity
## 19.2% Shared influence

############################################################################################
##                                                                                        ##
##                        Account for Spider Body Size and Country                        ##
##                                                                                        ## 
############################################################################################

varpart_mod4 <- varpart(dSpider,
                        ~ plot_data$Catchment + plot_data$Riparian + plot_data$Algae + plot_data$dis3 + plot_data$life2 + plot_data$wnb3 + plot_data$wnb5 + plot_data$SpiderCWM,
                        ~ Spider_FA_data_log$Body_Size,
                        ~ Spider_FA_data_log$Country)

plot(varpart_mod4)

## Include Country and Spider body size (Taxa level, not community)
## Account for Spider Identity - selected parameters
Gmm1 <- model.matrix( ~ Spider_body_size + dis3  + life2 + wnb3 + wnb5 + chl_a_m2_day + 
                        Riparian_PC1 + Catchment_PC1, Spider_FA_data_log) 
Gmm2 <- model.matrix( ~ Body_Size, Spider_FA_data_log)
Gmm3 <- model.matrix( ~ Country, Spider_FA_data_log) 
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 + Gmm2 + Gmm3)
RsquareAdj(rda_result)
vif.cca(rda_result)
## $adj.r.squared
## 0.3710135

## Predictors
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~ Gmm1 +  Condition(Gmm2) + Condition(Gmm3))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.02927923

## Spider body size
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~  Gmm2 +  Condition(Gmm1) + Condition(Gmm3))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.03628771

## Country
rda_result <- rda(Spider_FA_data_log[,c(12:48)] ~  Gmm3 +  Condition(Gmm1) + Condition(Gmm2))
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
## $adj.r.squared
## 0.1218503

## 0.3710135 - (0.02927923 + 0.03628771 + 0.1218503)
## 0.1835963

## 37.1% Combined influence of variation in PUFA 
## 2.9% Predictors independent of spider body size and country
## 3.6% Spider body size independent of environment and country
## 12.2% Country independent of environment and spider body size
## 18.4% Shared influence


#############################################################################################
##                                                                                         ##
##                        ADDITIONAL ANALYSES ALL SITES - TAXA                             ##        
##                                                                                         ##    
#############################################################################################

## Spider data represents functional and numerical responses to the environment
## Ground spiders are more prevalent in unbuffered sites, and this may influence the FA data available
## That is, ground spiders may be more likely to be found in the unbuffered sites
## Their larger body size and active hunting mode may influence their PUFA content

## head(Spider_FA_data_new)

## First check the number of spider taxa used
## See that Tetragnathidae and Lycosidae are the only taxa present at all sites
## Group taxa by site type
Spider_FA_community_sites <- Spider_FA_data_new %>%
  group_by(Country,Type_1) %>%
  count(Genus_family)

## Spread to wide format
Spider_FA_community_sites_wide <- Spider_FA_community_sites %>% spread(Genus_family,n)
## Turn NA to zeros
Spider_FA_community_sites_wide[is.na(Spider_FA_community_sites_wide)] <- 0
## Check output
head(Spider_FA_community_sites_wide)


## Need to test if ground spiders were more prevalent at unbuffered sites
## Group taxa by Mode of hunting
Spider_FA_community_site_Mode <- Spider_FA_data_new %>%
  group_by(Country,Type_1,Type_2,Site_block,Site_no,Site_name) %>%
  count(Mode)
## Spread to wide format
Spider_FA_community_site_mode_wide <- Spider_FA_community_site_Mode %>% spread(Mode,n)
## Turn NA to zeros
Spider_FA_community_site_mode_wide[is.na(Spider_FA_community_site_mode_wide)] <- 0


## Test if ground hunting spiders were more prevalent at unbuffered sites
## Notes about ANOVA https://stats.stackexchange.com/questions/60362/choice-between-type-i-type-ii-or-type-iii-anova
M1.1 <- glm(Ground ~ Type_1*Country, data=Spider_FA_community_site_mode_wide, family=poisson(link = "log"))
M1.2 <- glm(Ground ~ Type_1 + Country, data=Spider_FA_community_site_mode_wide, family=poisson(link = "log"))
BIC(M1.1)
BIC(M1.2)
sjPlot::tab_model(M1.2)
car::Anova(M1.2, "II")
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M1.2, ~ Type_1))
test(L.S, adjust = "tukey")
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M1.2, ~ Country))
test(L.S, adjust = "tukey")

## Test if web building spiders were more prevalent at buffered sites
M2.1 <- glm(Web ~ Type_1*Country, data=Spider_FA_community_site_mode_wide, family=poisson(link = "log"))
M2.2 <- glm(Web ~ Type_1 + Country, data=Spider_FA_community_site_mode_wide, family=poisson(link = "log"))
BIC(M2.1)
BIC(M2.2)
sjPlot::tab_model(M2.2)
car::Anova(M2.2, "II")
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M2.2, ~ Type_1))
test(L.S, adjust = "tukey")
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M2.2, ~ Country))
test(L.S, adjust = "tukey")

## Assess the mean counts for ground spiders as another indication of where differences occur
Spider_FA_community_site_mode_wide %>%
  group_by(Type_1) %>%
  summarise(mean=mean(Ground),
            sd=sd(Ground))

## Check if pisaurid spiders were more prevalent at unbuffered sites
## Group taxa
Spider_FA_community_sites <- Spider_FA_data_new %>%
  group_by(Country,Type_1,Type_2,Site_block,Site_no,Site_name) %>%
  count(Genus_family)
## Spread taxa to wide format
Spider_FA_community_sites_wide <- Spider_FA_community_sites %>% spread(Genus_family,n)
## Turn NA to zeros
Spider_FA_community_sites_wide[is.na(Spider_FA_community_sites_wide)] <- 0
## Check output
head(Spider_FA_community_sites_wide)

## Test if Pisauridae were more prevalent at unbuffered sites
## Model does not work with nested term for site_block
M3.1 <- glm(Pisauridae ~ Type_1*Country, data=Spider_FA_community_sites_wide, family=poisson(link = "log"))
M3.2 <- glm(Pisauridae ~ Type_1 + Country, data=Spider_FA_community_sites_wide, family=poisson(link = "log"))
BIC(M3.1)
BIC(M3.2)
sjPlot::tab_model(M3.2)
car::Anova(M3.1, "III") ## No interaction
car::Anova(M3.2, "II") ## Shows difference
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M3.2, ~ Type_1))
test(L.S, adjust = "tukey")
## Post-hoc contrasts to investigate differences
L.S <- pairs(lsmeans(M3.2, ~ Country))
test(L.S, adjust = "tukey")

## Assess the mean counts for ground spiders as another indication of where differences occur
Spider_FA_community_sites_wide %>%
  group_by(Type_1) %>%
  summarise(mean=mean(Pisauridae),
            sd=sd(Pisauridae))

#############################################################################################
##                               PERMANOVA ALL SITES - TAXA                                ##        
#############################################################################################

## Need to test how spiders used for FA analyses differed in terms of composition
## Use vegdist to create BC dissmilarity 
dspider <- decostand(Spider_FA_community_sites_wide[,c(8:13)], "hel")

## Test effect of site on composition of spiders used for FA after accounting for Country and Site blocks
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_community_sites_wide, Site_block)
adonis2(dspider ~ Country / Site_block + Type_1, method = "euc", by = "terms", data = Spider_FA_community_sites_wide, permutations = perm)
Mod0.1 <- adonis2(dspider ~ Country / Site_block + Type_1, method = "euc", by = "terms", data = Spider_FA_community_sites_wide, permutations = perm)

## Use vegdist to create Euclidean distance metric 
dspider_euc <- vegdist(dspider, method = "euc")

## Note that model cannot handle imbalance in site blocks (i.e., paired sites vs forest singles)
## Run Hellinger Euclidean pairwise PERMANOVA to test differences in site and country and interaction
permanova_0.1 <- pairwise.adonis2(dspider_euc ~ Type_1 + Country/Site_block, data = Spider_FA_community_sites_wide)
permanova_0.1

permanova_0.2 <- pairwise.adonis2(dspider_euc ~ Type_1 + Country, data = Spider_FA_community_sites_wide)
permanova_0.2

#############################################################################################
##                                                                                         ##
##                                 ALL SITES - FATTY ACIDS                                 ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                                     ALL SITES - TAXA                                    ##        
##                                                                                         ##    
#############################################################################################

## Create new data frame for transformations
Spider_FA_data_log <- Spider_FA_data_new

## Transform percentages using asin sqrt transformation
Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Transform body size using log transformation
Spider_FA_data_log[,c(84)] <- log(Spider_FA_data_log[,c(84)])

## BRAY-CURTIS

## Test interaction of taxa and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Taxa_code*Type_1, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.1 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Taxa_code*Type_1, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.1

## Test interaction of body size and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Body_Size, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.2 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Body_Size, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.2

## Test interaction of mode of hunting and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Mode + Taxa_code, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.3 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Mode + Taxa_code, method = "bray", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod1.3

## Tidyup workspace
rm(Mod1.1,Mod1.2,Mod1.3)

## EUCLIDEAN

## Test interaction of taxa and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Taxa_code*Type_1, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.1 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Taxa_code*Type_1, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.1

## Test interaction of body size and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Body_Size, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.2 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Body_Size, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.1

## Test interaction of mode of hunting and site
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)
adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Mode + Taxa_code, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.3 <- adonis2(Spider_FA_data_log[,c(13:49)] ~ Country / Site_block + Type_1 + Mode + Taxa_code, method = "euc", by = "terms", data = Spider_FA_data_log, permutations = perm)
Mod2.3

## Tidyup workspace
rm(Mod2.1,Mod2.2,Mod2.3)

#############################################################################################
##                                                                                         ##
##                                     SIMPER                                              ##        
##                                                                                         ##    
#############################################################################################

## Use SIMPER to test the importance of different FAs to changes among sites
## Note that SIMPER use Bray-Curtis dissimilarity to determine which species contribute to site differences
sim_all <- simper(Spider_FA_data_log[,c(13:49)], group=Spider_FA_data_log$Type_1, permutations = perm)

## Write summary object
sim_out <- summary(sim_all)

## Create summary object data frame
sim_out <- do.call(rbind.data.frame, sim_out)

## Sort species in order of decreasing significance to site differences
sim_out <- sim_out[order(sim_out$p),]

## Check output
head(sim_out[c(1:7),])

## Tidyup workspace
rm(sim_all,sim_out)

## Use SIMPER to test the importance of different FAs to changes among modes of hunting
## Note that SIMPER use Bray-Curtis dissimilarity to determine which species contribute to site differences
sim_all <- simper(Spider_FA_data_log[,c(13:49)], group=Spider_FA_data_log$Mode, permutations = 999)

## Write summary object
sim_out <- summary(sim_all)

## Create summary object data frame
sim_out <- do.call(rbind.data.frame, sim_out)

## Sort species in order of decreasing significance to site differences
sim_out <- sim_out[order(sim_out$p),]

## Check output
head(sim_out[c(1:7),])

## Tidyup workspace
rm(sim_all,sim_out)

#############################################################################################
##                                                                                         ##
##                        PAIRWISE ADONIS - TESTS EACH SITE CONTRAST                       ##        
##                                                                                         ##    
#############################################################################################

## BRAY-CURTIS

## Use vegdist to create BC dissimilarity 
dspider <- vegdist(Spider_FA_data_log[,c(13:49)], method = "bray")

head(Spider_FA_data_log)

## Run Bray-Curtis pairwise PERMANOVA to test interaction of taxa and site
permanova_1.1 <- pairwise.adonis2(dspider ~ Type_1*Taxa_code + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_1.1

## Run Bray-Curtis pairwise PERMANOVA to test interaction of size and site
permanova_1.2 <- pairwise.adonis2(dspider ~ Type_1*Body_Size + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_1.2

## Run Bray-Curtis pairwise PERMANOVA to test interaction of mode of hunting and site
permanova_1.3 <- pairwise.adonis2(dspider ~ Type_1*Mode + Taxa_code + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_1.3

## Tidyup workspace
rm(permanova_1.1,permanova_1.2,permanova_1.3)

## EUCLIDEAN

## Use vegdist to create Euclidean dissimilarity
dspider <- vegdist(Spider_FA_data_log[,c(13:49)], method = "euc")

## Run Euclidean pairwise PERMANOVA to test interaction of taxa and site
permanova_2.1 <- pairwise.adonis2(dspider ~ Type_1*Taxa_code + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_2.1

## Run Euclidean pairwise PERMANOVA to test interaction of size and site
permanova_2.2 <- pairwise.adonis2(dspider ~ Type_1*Body_Size + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_2.2

## Run Euclidean pairwise PERMANOVA to test interaction of mode of hunting and site
permanova_2.3 <- pairwise.adonis2(dspider ~ Type_1*Mode + Taxa_code + Country / Site_block, strata = 'Site_block', data = Spider_FA_data_log)
permanova_2.3

## Tidyup workspace
rm(permanova_2.1,permanova_2.2,permanova_2.3)

#############################################################################################
##                                                                                         ##
##                           BETADISPER - HOMOGENEITY OF DISPERSIONS                       ##        
##                                                                                         ##    
#############################################################################################

## Create new dataframe for plots
Spider_FA_data_log_plot <- Spider_FA_data_log

## Relable sites so that they appear in a desired order 
Spider_FA_data_log_plot[which(Spider_FA_data_log$Type_1=="FOR"),4]  <- "3_FOR"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Type_1=="UBF"),4]  <- "1_UBF"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Type_1=="FBF"),4]  <- "2_FBF"

## Relable sites so that they appear in a desired order 
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Agel"),3]  <- "5_Agel"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Aran"),3]  <- "4_Aran"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Liny"),3]  <- "6_Liny"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Lyco"),3]  <- "3_Lyco"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="OP"),3]    <- "1_OP"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Pisa"),3]  <- "2_Pisa"
Spider_FA_data_log_plot[which(Spider_FA_data_log$Taxa_code=="Tetr"),3]  <- "7_Tetr"

## Create new dataframe for plots
Spider_FA_data_log_plot2 <- Spider_FA_data_log

## Relable sites so that they appear in a desired order 
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Agelenidae"),3]  <- "Ageli"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Araneidae"),3]  <- "Arane"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Linyphiidae"),3]  <- "Linyp"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Lycosidae"),3]  <- "Lycos"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Opiliones"),3]    <- "Opili"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Pisauridae"),3]  <- "Pisau"
Spider_FA_data_log_plot2[which(Spider_FA_data_log_plot$Genus_family=="Tetragnathidae"),3]  <- "Tetra"

## BRAY-CURTIS

## Use vegdist to create BC dissmilarity 
dspider <- vegdist(Spider_FA_data_log_plot[,c(13:49)], method = "bray")

## Plot feeding traits and other key traits and Gammarid abundances 
png(filename="Burdon_et_al_Fig.S1.png", 
    type="cairo",
    units="in", 
    width=11, 
    height=13, 
    pointsize=16, 
    res=600)

## Countries, Taxa, Mode, Site
par(mfrow=c(4,2), mar = c(4,4,2,2) + 0.1)

## Set permutations
perm <- how(nperm = 999)
setBlocks(perm) <- with(Spider_FA_data_log, Site_block)



## Differences in dispersion of: Site type
beta_FA <- betadisper(dspider, Spider_FA_data_log$Type_1)
str(beta_FA)
plot(beta_FA, main="Site", col =c("Green","Black","Red"), pch=c(1,2,3))
legend("topleft", legend=c("FBF","FOR","UBF"), pch = c(1,2,3), col =c("Green","Black","Red"))

beta_FA <- betadisper(dspider, Spider_FA_data_log_plot$Type_1)
boxplot(beta_FA, xlab="Site", xaxt="n")
axis(1, at=1:3, labels=c("Unbuffered","Buffered","Forest"))

# F-test
anova(beta_FA)
TukeyHSD(beta_FA)
betadisper_posthoc_bray_site <- TukeyHSD(beta_FA)

# Permutation test
permutest(beta_FA, permutations = perm)
betadisper_permute_bray_site <- permutest(beta_FA, permutations = perm)



## Differences in dispersion of: Mode of hunting
beta_FA <- betadisper(dspider, Spider_FA_data_log_plot$Mode)
plot(beta_FA, main="Mode of hunting", col =c("Black","firebrick4"), pch=c(1,2))
legend("topleft", legend=c("Ground","web"), pch = c(1,2), col =c("Black","firebrick4"))

boxplot(beta_FA, xlab="Mode of hunting") # plot how plots are spread within locations

# F-test
anova(beta_FA)
TukeyHSD(beta_FA)
betadisper_posthoc_bray_mode <- TukeyHSD(beta_FA)

# Permutation test
permutest(beta_FA, permutations = perm)
betadisper_permute_bray_mode <- permutest(beta_FA, permutations = perm)



## Differences in dispersion of: Spider taxa
beta_FA <- betadisper(dspider, Spider_FA_data_log_plot2$Taxa_code)
plot(beta_FA, main="Taxa", col =c("Black","lightsalmon","springgreen2","steelblue1","cyan1","violetred","gold1"), pch=c(1:7))
legend("topleft", legend=c("Ageli","Arane","Linyp","Lycos","Opili","Pisau","Tetra"), pch = c(1:7), col =c("Black","lightsalmon","springgreen2","steelblue1","cyan1","violetred","gold1"))


beta_FA <- betadisper(dspider, Spider_FA_data_log_plot$Taxa_code)
boxplot(beta_FA, xlab="Taxa", xaxt="n") # plot how plots are spread within locations
axis(1, at=1:7, labels=c("Opili","Pisau","Lycos","Arane",
                         "Ageli","Linyp","Tetra"))

# F-test
anova(beta_FA)
TukeyHSD(beta_FA)
betadisper_posthoc_bray_taxa <- TukeyHSD(beta_FA)

# Permutation test
permutest(beta_FA, permutations = perm)
betadisper_permute_bray_taxa <- permutest(beta_FA, permutations = perm)



## Differences in dispersion of: Countries
beta_FA <- betadisper(dspider, Spider_FA_data_log_plot$Country)
plot(beta_FA, main="Country", col =c("Black","maroon","springgreen4","steelblue3"), pch=c(1:7))
legend("topleft", legend=c("BE","NO","RO","SE"), pch = c(1:7), col =c("Black","maroon","springgreen4","steelblue3"))

boxplot(beta_FA, xlab="Country") # plot how plots are spread within locations

# F-test
anova(beta_FA)
TukeyHSD(beta_FA)
betadisper_posthoc_bray_country <- TukeyHSD(beta_FA)

# Permutation test
permutest(beta_FA, permutations = perm)
betadisper_permute_bray_country <- permutest(beta_FA, permutations = perm)

## Close plot
dev.off()

## Tidyup workspace
rm(betadisper_posthoc_bray_country, betadisper_permute_bray_country,
   betadisper_posthoc_bray_taxa, betadisper_permute_bray_taxa,
   betadisper_posthoc_bray_mode, betadisper_permute_bray_mode,
   betadisper_posthoc_bray_site, betadisper_permute_bray_site)


#############################################################################################
##                                                                                         ##
##                                 SPIDERS BY GROUPS                                       ## 
##                            1. Mode of hunting (ground, web)                             ##
##                            2. Key families (Lycosidae, Tetragnathidae)                  ##
##                                                                                         ##    
#############################################################################################

## Create data frames for different modes of hunting
Spider_FA_data_ground <- Spider_FA_data_new[which(Spider_FA_data_new$Mode=="Ground"),]
Spider_FA_data_web <- Spider_FA_data_new[which(Spider_FA_data_new$Mode=="Web"),]

## Create data frames for key spider Families
Spider_FA_data_lycos <- Spider_FA_data_new[which(Spider_FA_data_new$Genus_family=="Lycosidae"),]
Spider_FA_data_tetra <- Spider_FA_data_new[which(Spider_FA_data_new$Genus_family=="Tetragnathidae"),]

#############################################################################################
##                                                                                         ##
##                                 GROUND HUNTING SPIDERS                                  ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset
Spider_FA_data_log <- Spider_FA_data_ground

## Transform percentages using asin sqrt transformation
Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- na.omit(Spider_FA_data_log[,-12])

## Remove NAs (256 -> 240 samples)
Spider_FA_data_log <- na.omit(Spider_FA_data_log)

head(Spider_FA_data_log)

## Check the number of sites used
Data_sites <- Spider_FA_data_log %>%
  group_by(Country)%>%
  count(Site_name)

## Create data frame and add one for each site
Data_sites$Site <- 1

## Count sites from each Country
Data_sites %>%
  group_by(Country)%>%
  count(Site)

colnames(Spider_FA_data_log)

## Transformations
Spider_FA_data_log[,c(50)] <- log1p(Spider_FA_data_log[,c(50)]) ## Log+1 transform chl-a concentrations per day
Spider_FA_data_log[,c(51)] <- log(Spider_FA_data_log[,c(51)]) ## Log transform stream wetted width (m)
Spider_FA_data_log[,c(87)] <- log(Spider_FA_data_log[,c(87)]) ## Log-transform body size

## Standardization
Spider_FA_data_log[,c(49:87)] <- decostand(Spider_FA_data_log[,c(49:87)],"standardize",na.action=na.exclude)

## Create plot data
colnames(Spider_FA_data_log)
plot_data <- Spider_FA_data_log[,c(50,52,61,71,84,86,87)]
colnames(plot_data)
names(plot_data) <- c("Algae","Riparian","cd1","life2","wnb3","wnb5","SpiderBS") 

## Create data frame for PUFA data
dSpider <- Spider_FA_data_log[,c(12:48)]

FA_lookup$Specific_PUFA[c(13,14,17,18,24)] 

names(dSpider)[13] <- "ALA"
#names(dSpider)[14] <- "LIN"
#names(dSpider)[17] <- "ARA"
names(dSpider)[18] <- "EPA"
names(dSpider)[24] <- "DHA"

# Analyse and plot the RDA using ggplot (ggord package)
dSpider <- Spider_FA_data_log[,c(12:48)]
colnames(dSpider) <- c(FA_lookup$Short_name_2[-38])

rda_result <- rda(dSpider ~ SpiderBS + Algae + Riparian + cd1 + life2 + wnb3 + wnb5, plot_data) 
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
anova.cca(rda_result, by="axis")
vif.cca(rda_result)

p1 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "2"), ellipse = TRUE,
      ellipse_pro = 0.95, ptslab=T, obslab=F, size= 6, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    guide=F) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), 
                      guide=F) +
  guides(fill=guide_legend(title="Country")) + 
  #expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

#############################################################################################
##                                                                                         ##
##                                  WEB-BUILDING SPIDERS                                   ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset
Spider_FA_data_log <- Spider_FA_data_web

## Transform percentages using asin sqrt transformation
Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- na.omit(Spider_FA_data_log[,-12])

## Remove NAs (256 -> 240 samples)
Spider_FA_data_log <- na.omit(Spider_FA_data_log)

head(Spider_FA_data_log)

## Check the number of sites used
Data_sites <- Spider_FA_data_log %>%
  group_by(Country)%>%
  count(Site_name)

## Create data frame and add one for each site
Data_sites$Site <- 1

## Count sites from each Country
Data_sites %>%
  group_by(Country)%>%
  count(Site)

colnames(Spider_FA_data_log)

## Transformations
Spider_FA_data_log[,c(50)] <- log1p(Spider_FA_data_log[,c(50)]) ## Log+1 transform chl-a concentrations per day
Spider_FA_data_log[,c(51)] <- log(Spider_FA_data_log[,c(51)]) ## Log transform stream wetted width (m)
Spider_FA_data_log[,c(87)] <- log(Spider_FA_data_log[,c(87)]) ## Log-transform body size

## Standardization
Spider_FA_data_log[,c(49:87)] <- decostand(Spider_FA_data_log[,c(49:87)],"standardize",na.action=na.exclude)

## Check output
colnames(Spider_FA_data_log)
plot_data <- Spider_FA_data_log[,c(49,50,53,71,84,86,87)]
colnames(plot_data)

names(plot_data) <- c("SpiderCWM","Algae","Catchment","life2","wnb3","wnb5","SpiderBS") 

## Create data frame for PUFA data
dSpider <- Spider_FA_data_log[,c(12:48)]

head(dSpider)

FA_lookup$Specific_PUFA[c(13,14,17,18,24)] 

names(dSpider)[13] <- "ALA"
#names(dSpider)[14] <- "LIN"
#names(dSpider)[17] <- "ARA"
names(dSpider)[18] <- "EPA"
names(dSpider)[24] <- "DHA"

# Analyse and plot the RDA using ggplot (ggord package)
dSpider <- Spider_FA_data_log[,c(12:48)]
colnames(dSpider) <- c(FA_lookup$Short_name_2[-38])

rda_result <- rda(dSpider ~ SpiderCWM + Algae + Catchment + life2 + wnb3 + wnb5 + SpiderBS, plot_data) 
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
anova.cca(rda_result, by="axis")
vif.cca(rda_result)

p2 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "2"), ellipse = TRUE,
      ellipse_pro = 0.95, ptslab=T, obslab=F, size= 6, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    guide=F) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), 
                      guide=F) +
  guides(fill=guide_legend(title="Country")) + 
  #expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

#############################################################################################
##                                                                                         ##
##                                      LYCOSIDAE                                          ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset
Spider_FA_data_log <- Spider_FA_data_lycos 

## Transform percentages using asin sqrt transformation
Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- na.omit(Spider_FA_data_log[,-12])

## Remove NAs (256 -> 240 samples)
Spider_FA_data_log <- na.omit(Spider_FA_data_log)

head(Spider_FA_data_log)

## Check the number of sites used
Data_sites <- Spider_FA_data_log %>%
  group_by(Country)%>%
  count(Site_name)

## Create data frame and add one for each site
Data_sites$Site <- 1

## Count sites from each Country
Data_sites %>%
  group_by(Country)%>%
  count(Site)

colnames(Spider_FA_data_log)

## Transformations
Spider_FA_data_log[,c(50)] <- log1p(Spider_FA_data_log[,c(50)]) ## Log+1 transform chl-a concentrations per day
Spider_FA_data_log[,c(51)] <- log(Spider_FA_data_log[,c(51)]) ## Log transform stream wetted width (m)
Spider_FA_data_log[,c(87)] <- log(Spider_FA_data_log[,c(87)]) ## Log-transform body size

## Standardization
Spider_FA_data_log[,c(49:87)] <- decostand(Spider_FA_data_log[,c(49:87)],"standardize",na.action=na.exclude)

plot_data <- Spider_FA_data_log[,c(50,52,61,71,84,86)]
colnames(plot_data)

names(plot_data) <- c("Algae","Riparian","cd1","life2","wnb3","wnb5") 

## Create data frame for PUFA data
#dSpider <- Spider_FA_data_log[,c(12:48)]
#head(dSpider)

#FA_lookup$Specific_PUFA[c(13,14,17,18,24)] 

#names(dSpider)[13] <- "ALA"
#names(dSpider)[14] <- "LIN"
#names(dSpider)[17] <- "ARA"
#names(dSpider)[18] <- "EPA"
#names(dSpider)[24] <- "DHA"

# Analyse and plot the RDA using ggplot (ggord package)
dSpider <- Spider_FA_data_log[,c(12:48)]
colnames(dSpider) <- c(FA_lookup$Short_name_2[-38])

rda_result <- rda(dSpider ~ Algae + Riparian + cd1 + life2 + wnb3 + wnb5, plot_data) 
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
anova.cca(rda_result, by="axis")
vif.cca(rda_result)

require(ggord)

p3 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "2"), ellipse = TRUE,
      ellipse_pro = 0.95, ptslab=T, obslab=F, size= 6, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    guide=F) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), 
                      guide=F) +
  guides(fill=guide_legend(title="Country")) + 
  #expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

#############################################################################################
##                                                                                         ##
##                                    TETRAGNATHIDAE                                       ##        
##                                                                                         ##    
#############################################################################################

## Create dummy dataset
Spider_FA_data_log <- Spider_FA_data_tetra

## Transform percentages using asin sqrt transformation
Spider_FA_data_log[,c(13:49)] <- asin(sqrt(Spider_FA_data_log[,c(13:49)]/100))

## Remove Tribe (not applicable to spiders)
Spider_FA_data_log <- na.omit(Spider_FA_data_log[,-12])

## Remove NAs (256 -> 240 samples)
Spider_FA_data_log <- na.omit(Spider_FA_data_log)

head(Spider_FA_data_log)

## Check the number of sites used
Data_sites <- Spider_FA_data_log %>%
  group_by(Country)%>%
  count(Site_name)

## Create data frame and add one for each site
Data_sites$Site <- 1

## Count sites from each Country
Data_sites %>%
  group_by(Country)%>%
  count(Site)

colnames(Spider_FA_data_log)

## Transformations
Spider_FA_data_log[,c(50)] <- log1p(Spider_FA_data_log[,c(50)]) ## Log+1 transform chl-a concentrations per day
Spider_FA_data_log[,c(51)] <- log(Spider_FA_data_log[,c(51)]) ## Log transform stream wetted width (m)
Spider_FA_data_log[,c(87)] <- log(Spider_FA_data_log[,c(87)]) ## Log-transform body size

## Standardization
Spider_FA_data_log[,c(49:87)] <- decostand(Spider_FA_data_log[,c(49:87)],"standardize",na.action=na.exclude)

## Check output
colnames(Spider_FA_data_log)

## Create plot data for selected drivers
plot_data <- Spider_FA_data_log[,c(52,53,84,86)]
head(plot_data)

names(plot_data) <- c("Riparian","Catchment","wnb3","wnb5") 


## Create data frame for PUFA data
#dSpider <- Spider_FA_data_log[,c(12:48)]

head(dSpider)

#FA_lookup$Specific_PUFA[c(13,14,17,18,24)] 

#names(dSpider)[13] <- "ALA"
#names(dSpider)[14] <- "LIN"
#names(dSpider)[17] <- "ARA"
#names(dSpider)[18] <- "EPA"
#names(dSpider)[24] <- "DHA"

# Analyse and plot the RDA using ggplot (ggord package)
dSpider <- Spider_FA_data_log[,c(12:48)]
colnames(dSpider) <- c(FA_lookup$Short_name_2[-38])

rda_result <- rda(dSpider ~ Riparian + Catchment +  wnb3 + wnb5, plot_data)  
RsquareAdj(rda_result)
anova(rda_result, step=999, perm.max=999)
anova.cca(rda_result, by="axis")
vif.cca(rda_result)

p4 <- ggord(rda_result, Spider_FA_data_log$Country, axes = c("1", "2"), ellipse = TRUE,
      ellipse_pro = 0.95, ptslab=T, obslab=F, size= 6, alpha=0.6, repel = TRUE) +
  scale_fill_manual(breaks=c("RO", "SE", "NO", "BE"),
                    values=c("forestgreen","olivedrab4","tan4","grey50"),
                    labels = c("Romania", "Sweden", "Norway", "Belgium"),
                    guide=F) +
  scale_colour_manual(breaks=c("RO", "SE", "NO", "BE"),
                      values=c("forestgreen","olivedrab4","tan4","grey50"),
                      labels = c("Romania", "Sweden", "Norway", "Belgium"), 
                      guide=F) +
  guides(fill=guide_legend(title="Country")) + 
  #expand_limits(x=c(-0.8,1.2)) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 15, face="bold"),
        axis.title.x=element_text(vjust=-0.65, size = 15, face="bold"),
        axis.text.y=element_text(size = 13),
        axis.text.x=element_text(size = 13))

#############################################################################################
##                                                                                         ##
##                                  CREATE COMBINED PLOT                                   ##        
##                                                                                         ##    
#############################################################################################

png(filename="Burdon_et_al_Fig.S3.png", 
    type="cairo",
    units="in", 
    width=9, 
    height=8, 
    pointsize=20, 
    res=600)

## Arrange using grid plot
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="right")
grid.text("(a)", x = unit(0.03, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(c)", x = unit(0.03, "npc"), y = unit(0.47, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(b)", x = unit(0.47, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(d)", x = unit(0.47, "npc"), y = unit(0.47, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))

dev.off()

#############################################################################################
##                                                                                         ##
##                                VARIATION PARTITIONING                                   ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                                 CREATE COMBINED PLOT                                    ##        
##                                                                                         ##    
#############################################################################################

png(filename="Burdon_et_al_Fig.S4.png", 
    type="cairo",
    units="in", 
    width=16, 
    height=10, 
    pointsize=20, 
    res=600)

par(mfrow=c(2,4)) 

## GROUND
## Ground: Main hypothesis
MyVenn <- venneuler(c(Spider=0.04877, Insect=0.11720, Algae=0.02922, Environment=0.09219,
                      "Spider&Insect"=0.00969,
                      "Insect&Algae"=0.00015,
                      "Spider&Algae"=0,
                      "Spider&Environment"=0,
                      "Environment&Insect"=0,
                      "Environment&Algae"=0,
                      "Spider&Insect&Environment"=0,
                      "Spider&Insect&Algae"=0.00496,
                      "Environment&Insect&Algae"=0.00197,
                      "Spider&Algae&Environment"=0.00112,
                      "Spider&Insect&Algae&Environment"=0.00040 
))

MyVenn$labels <- c("","","","")

plot(MyVenn, col=c('grey50','peru','yellowgreen','olivedrab'))

## Ground: H2
MyVenn <- venneuler(c(Riparian=0.06480, Foodweb=0.15482, Environment=0,
                      "Riparian&Foodweb"=0.00219,
                      "Foodweb&Environment"=0.05194,
                      "Riparian&Environment"=0.02739,
                      "Riparian&Foodweb&Environment"=0))

MyVenn$labels <- c("","","")

plot(MyVenn, col=c('forestgreen','peru','olivedrab'))

## WEB
## Web: Main hypothesis
MyVenn <- venneuler(c(Spider=0.03916, Insect=0.03895, Algae=0.01159, Environment=0.01745,
                      "Spider&Insect"=0.00369,
                      "Insect&Algae"=0.00230,
                      "Spider&Algae"=0.00193,
                      "Spider&Environment"=0.00518,
                      "Environment&Insect"=0.05651,
                      "Environment&Algae"=0,
                      "Spider&Insect&Environment"=0.01029,
                      "Spider&Insect&Algae"=0,
                      "Environment&Insect&Algae"=0.02174,
                      "Spider&Algae&Environment"=0.00194,
                      "Spider&Insect&Algae&Environment"=0.01701 
))

MyVenn$labels <- c("","","","")

plot(MyVenn, col=c('grey50','peru','yellowgreen','olivedrab'))

## Web: H2
MyVenn <- venneuler(c(Riparian=0.00330, Foodweb=0.08680, Environment=0.01398,
                      "Riparian&Foodweb"=0.00915,
                      "Foodweb&Environment"=0.04152,
                      "Riparian&Environment"=0.00348,
                      "Riparian&Foodweb&Environment"=0.07014))

MyVenn$labels <- c("","","")

plot(MyVenn, col=c('forestgreen','peru','olivedrab'))

## LYCOSIDAE
## Lycosidae: Main hypothesis
MyVenn <- venneuler(c(Insect=0.15514, Algae=0.02046, Environment=0.14598,
                      "Insect&Algae"=0.00981,
                      "Environment&Algae"=0,
                      "Environment&Insect"=0,
                      "Insect&Algae&Environment"=0.00358 
))

MyVenn$labels <- c("","","","")

plot(MyVenn, col=c('peru','yellowgreen','olivedrab'))

## Lycosidae: H2
MyVenn <- venneuler(c(Riparian=0.10167, Foodweb=0.13201, Environment=0,
                      "Riparian&Foodweb"=0.01135,
                      "Foodweb&Environment"=0.05340,
                      "Riparian&Environment"=0.04431,
                      "Riparian&Foodweb&Environment"=0))

MyVenn$labels <- c("","","")

plot(MyVenn, col=c('forestgreen','peru','olivedrab'))

## TETRAGNATHIDAE
## Tetragnathidae: Main hypothesis
MyVenn <- venneuler(c(Insect=0.05300, Environment=0.06848,
                      "Insect&Environment"=0.10511
))

MyVenn$labels <- c("","")

plot(MyVenn, col=c('peru','olivedrab'))

## Tetragnathidae: H2
MyVenn <- venneuler(c(Riparian=0.02360, Foodweb=0.05300, Environment=0.04601,
                      "Riparian&Foodweb"=0,
                      "Foodweb&Environment"=0.04863,
                      "Riparian&Environment"=0.03549,
                      "Riparian&Foodweb&Environment"=0.02320))

MyVenn$labels <- c("","","")

plot(MyVenn, col=c('forestgreen','peru','olivedrab'))

dev.off()
