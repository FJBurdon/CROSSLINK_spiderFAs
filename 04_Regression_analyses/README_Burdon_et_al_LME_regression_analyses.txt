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
## 15.    s3: Maximum body size (cm) ≥ 0.5–1
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