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