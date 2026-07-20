###############################################################################################################
##                                                                                                           ##
##   DATA ANALYSIS relating to Fig.6 in "Fatty acid biomarkers reveal landscape influences on               ##
##   linkages between aquatic and terrestrial food webs" by Burdon et al. (in press) Ecological Monographs   ##        
##                                                                                                           ##    
###############################################################################################################

# Analysis of CROSSLINK data from field survey (All case-study basins)
# Version: 1.0
# Author: Dr. Francis J. Burdon
# Started: 10 July 2019 
# Updated: 7 June 2025 

## This R script performs the correlation analysis of 10 community-weighted mean (CWM) invertebrate dispersal traits
## This trait data is based on site stream macroinvertebrate data from six pooled Surber samples (3 erosional, 3 depositional) 
## These 10 aquatic invertebrate dispersal traits from Sarremejane et al. (2020) were hypothesized to be important for 
## aquatic-terrestrial trophic connectivity in temperate streams.

## Sarremejane, R., et al. (2020). DISPERSE, a trait database to assess the dispersal potential of European aquatic macroinvertebrates. 
## Scientific Data 7:386. doi:10.1038/s41597-020-00732-7


## Metadata for "Burdon_et_al_Fig.6_invert_dispersal_trait_data.csv" (Numbers = column)
## 1.     Site: Site number	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Site_identifier: Unique code given to each site based on country, site, and site type 
## 4.     Site_code_2: Unique code given to each site in each case-study basin
## 5.     Site_name: Site name	
## 6.     Type_1:	Site type (forest, unbuffered, buffered)
## 7.     Type_2:	Site type (paired or reference)
## 8.     Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 9.     s3: Maximum body size (cm) ≥ 0.5–1
## 10.    cd1: Life cycle duration ≤ 1 year
## 11.    cy2: Potential number of reproductive cycles per year - One
## 12.    dis3: Dispersal strategy Aerial passive
## 13.    dis4: Dispersal strategy Aerial active
## 14.    life1: Adult life span < 1 week
## 15.    life2: Adult life span ≥ 1 week – 1 month
## 16.    wnb2: Wing pair type - 1 pair + halters
## 17.    wnb3: Wing pair type - 1 pair + 1 pair of small hind wings
## 18.    wnb5: Wing pair type - 2 similar-sized pairs