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