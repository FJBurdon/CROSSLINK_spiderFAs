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