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

## This scripts performs analyses and creates figures relating to multivariate analyses of spider fatty acid composition
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
## 54.    Catchment_PC1 : Catchment impact PC1 (42.3%)	
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