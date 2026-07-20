#########################################################################################################################
##                                                                                                                     ##  
##                                    README: DATA PACKAGE FOR BURDON ET AL.                                           ##
##  	Fatty acid biomarkers reveal landscape influences on linkages between aquatic and terrestrial food webs        ##
##                                            Ecological Monographs                                                    ##         
##                                                                                                                     ##  
#########################################################################################################################

Dataset Version and Release History
-----------------------------------

* File name: README_Fatty_Acid_Riparian_Foodwebs_v.1.0.txt
* Authors: Francis J. Burdon
* Other contributors: Jasmina Sargac, Ellinor Ramberg, Cristina Popescu, Nita Darmina, Corina Bradu, Marie A. E. Forio, Felix Witing, Benjamin Kupilas, Danny C. P. Lau, Martin Volk, Geta Rîşnoveanu, Peter Goethals, Nikolai Friberg, Richard K. Johnson, Brendan G. McKie  
* Date created: 2025-06-09
* Date modified: 2025-06-09

Dataset Version and Release History
-----------------------------------

* Current Version:
* Number: 1.1
* Date: 2025-06-12
* Persistent identifier: https://doi.org/10.5061/dryad.931zcrjxz
* Summary of changes: Updated title, added URL for dataset

* Embargo Provenance: n/a
  * Scope of embargo: n/a
  * Embargo period: n/a

Dataset Attribution and Usage
-----------------------------

* Dataset Title: Data for the article "Fatty acid biomarkers reveal landscape influences on linkages between aquatic and terrestrial food webs"

* Persistent Identifier: 

* Dataset Contributors:

  * Creators: Francis J. Burdon, Jasmina Sargac, Ellinor Ramberg, Cristina Popescu, Nita Darmina, Corina Bradu, Marie A. E. Forio, Felix Witing, Benjamin Kupilas, Danny C. P. Lau, Martin Volk, Geta Rîşnoveanu, Peter Goethals, Nikolai Friberg, Richard K. Johnson, Brendan G. McKie  

* Date of Issue: 2025-06-09

* License: Use of these data is covered by the following license:
  * Title: CC0 1.0 Universal (CC0 1.0)
  * Specification: https://creativecommons.org/publicdomain/zero/1.0/; the authors respectfully request to be contacted by researchers interested in the re-use of these data so that the possibility of collaboration can be discussed. 

* Suggested Citations:

  * Dataset citation:
    > Burdon, F. J., Sargac, J., Ramberg, E., Cristina, P., Darmina, N., Bradu, C., Forio, M. A. E., Witing, F., Kupilas, B., Lau, D. C. P., Volk, M., Rîşnoveanu, G., Goethals, P., Friberg, N., Johnson, R. K., McKie, B. G.   
 (2025). Data from: Fatty acid biomarkers reveal landscape influences on linkages between aquatic and terrestrial food webs. Dryad Digital Repository. https://doi.org/

  * Corresponding publication:
    > Burdon, F. J., Sargac, J., Ramberg, E., Cristina, P., Darmina, N., Bradu, C., Forio, M. A. E., Witing, F., Kupilas, B., Lau, D. C. P., Volk, M., Rîşnoveanu, G., Goethals, P., Friberg, N., Johnson, R. K., McKie, B. G. (2025). Fatty acid biomarkers reveal landscape influences on linkages between aquatic and terrestrial food webs. Ecological Monographs. Accepted. DOI: 

Contact Information
-------------------

  * Name: Francis J. Burdon
  * Affiliations: SLU, Sweden and School of Science, University of Waikato, Hamilton, NZ
  * ORCID ID: https://orcid.org/0000-0002-5398-4993
  * Email: francis.burdon@waikato.ac.nz
  * Alternate Email: f.burdon@gmail.com
  * Address: e-mail preferred

- - -

Additional Dataset Metadata
===========================

Acknowledgements
----------------

* Funding sources: This research for the CROSSLINK project was funded through the 2015-2016 BiodivERsA COFUND scheme. National funders were: the Swedish Research Council for Sustainable Development (FORMAS, project 2016-01945) and the Swedish Environmental Protection Agency; The Research Council of Norway (NFR, project 264499); The Research Foundation of Flanders (FWO Belgium, project G0H6516N); the Romanian National Authority for Scientific Research and Innovation (CCCDI – UEFISCDI, project BiodivERsA3-2015-49-CROSSLINK within PNCDI III); and the German Federal Ministry of Education and Research (BMBF, project FKZ: 01LC1621A). F.J.B. acknowledges funding from the Marsden Fund Council managed by Royal Society Te Apārangi (MFP-UOW2303).


Dates and Locations
-------------------

* Dates of data collection: Aquatic and terrestrial samples were collected between Autumn 2017 and Summer 2018 (see publication for more details)

* Geographic locations of data collection: Fieldwork conducted in Sweden, Norway, Romania, and Belgium (see publication for more details)

* Other locations pertaining to dataset contents: The fatty acid composition of riparian spiders was analysed at the Swedish Metabolomics Centre, Umeå, Sweden.

- - -

Methodological Information
==========================

* Methods of data collection/generation: see manuscript for details

- - -

Data and Folder Overview
======================

Summary Metrics
---------------

* Main Folder count: 6
* See below for specific information for data and scripts in each folder
* Data package size: 1.36 MB
* Script file formats: .R
* Data file formats: .csv

Table of Contents (Overview)
----------------------------

* NOTES: This data package contains six main folders:

1. "01_PCA_gradients" - Data from riparian and catchment PCAs and R scripts used to generate data used in subsequent analyses
2. "02_Multivariate_analyses" - Data and R scripts for testing determinants of fatty acid composition in riparian spiders using redundancy analysis and variation partitioning
3. "03_Dispersal_traits" - Data and R scripts for correlating selected aquatic invertebrate dispersal traits and providing evidence for trait syndromes
4. "04_Regression_analyses" - Data and R scripts for linear mixed effect (LME) models testing the correlation of realized trophic connectivity in riparian spiders with aquatic prey dispersal traits 
5. "05_SEM" - Data and R scripts for the structural equation models (SEMs) testing predicted direct and indirect effects on realized trophic connectivity   
6. "06_Algae_spider_FAs" - Data and R scripts for the generalized additive models (GAMs) used to investigate a 'subsidy-stress' relationship between stream periphyton and realized trophic connectivity in riparian spiders

Setup
-----

* Unpacking instructions: n/a

* Recommended software/tools: RStudio 2024.09.1; R version 4.4.2

- - -

More specific information for each of the above folders:

First folder: "01_PCA_gradients"
---------------------------------------

* NOTES: 
	* Data and R scripts used to generate Figure 3 and perform analyses testing differences between Site Type, Country and their interaction

Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 3
* Total file size: 23.6 KB
* Range of individual file sizes: 2 KB - 12 KB
* File formats: csv, .R, and .txt

Table of Contents
-----------------

RScript_catchment_riparian_gradients_analysis.R
Burdon_et_al_Fig.3_catchment_riparian_data.csv
README_Burdon_et_al_PCA_gradients.txt

Details for: "Burdon_et_al_Fig.3_catchment_riparian_data.csv"
---------------------------------------

* Description: A comma-delimited file containing principle component analysis (PCA) site scores for catchment and riparian attributes of stream sites in Sweden, Norway, Romania, and Belgium

* Format(s): .csv

* Size(s): 11 KB

* Dimensions: 103 rows by 9 columns

* Variables:
## 1.   Site_number: Site number at the country level
## 2.   Country: Country where data was collected (Sweden, Norway, Romania, Belgium)	
## 3.   Type_1:	Site type (forest, unbuffered, buffered)
## 4.   Type_2:	Site type (paired or reference)
## 5.   Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 6.   Catchment_PC1: Catchment impact PC1 (42.3%)	
## 7.   Catchment_PC2: Catchment impact PC2 (17.5%)	
## 8.   Riparian_PC1:	Riparian condition PC1 (36.7%)
## 9.   Riparian_PC2:	Riparian condition PC2 (14.2%)

* Missing data codes: NA


Second folder: "02_Multivariate_analyses"
---------------------------------------

* NOTES: 
* Includes input data and R script for testing determinants of fatty acid composition in riparian spiders using redundancy analysis and variation partitioning

Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 4
* Total file size: 456 KB
* Range of individual file sizes: 3 KB - 363 KB
* File formats: .csv, .R, and .txt

Table of Contents
-----------------

"RScript_multivariate_analyses_RDA_varpart.R"
"Burdon_et_al_Figs.4_5_multivariate_data.csv"
"Burdon_et_al_Figs.4_5_FA_lookup_table.csv"
"README_Burdon_et_al_multivariate_analyses.txt"


Details for:"Burdon_et_al_Figs.4_5_multivariate_data.csv"
---------------------------------------

* Description: A comma-delimited file containing multivariate data from stream sites in Sweden, Norway, Romania, and Belgium
* Includes 33 freshwater invertebrate dispersal traits from Sarremejane et al. (2020) 
* Sarremejane, R., et al. (2020). DISPERSE, a trait database to assess the dispersal potential of European aquatic macroinvertebrates. Scientific Data 7:386. doi:10.1038/s41597-020-00732-7

* Format(s): .csv

* Size(s): 363 KB

* Dimensions: 257 rows by 88 columns

* Variables:
## 1.     Genus_family: Family or taxonomic groups used for arachnids	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Taxa_code: Short name for arachnid family or taxonomic group	
## 4.     Type_1: Site type (forest, unbuffered, buffered)
## 5.     Type_2: Site type (paired or reference)
## 6.     Site_block: Blocking term for paired sites (unbuffered, buffered)
## 7.     Site_no:  Site number	
## 8.     Site_name:  Site name	
## 9.     Group:  Arachnida	
## 10.    Mode:  Mode of hunting for arachnids
## 11.    Lowest_level:  Lowest taxonomic level used	
## 12.    Tribe: Used for ground beetles not included in this dataset
## 13-49. FA1-FA37: Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
## 50.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 51.    chl_a_m2_day: Algal accrual on tiles
## 52.    Width: Stream wetted channel width (m)
## 53.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 54.    Catchment_PC1 :  Catchment impact PC1 (42.3%)	
## 55-87.  s1-wnb5: CWM abundances for 33 invertebrate dispersal trait modalities (e.g., aerial active, body size)
## 88.   Body_Size: Mean body size of spider Family sampled for FAs

* Missing data codes: NA


Details for:"Burdon_et_al_Figs.4_5_FA_lookup_table.csv"
---------------------------------------

* Description: A comma-delimited file containing fatty acid names and grouping information 

* Format(s): .csv

* Size(s): 3 KB

* Dimensions: 39 rows by 7 columns

* Variables:
## 1.     FA: Code for fatty acids (FA1-FA38)
## 2.     Short_name_2: Common chemical nomenclature used for fatty acids (FAs)
## 3.     FA_1: FA grouping 1 (separates out omega-3 and omega-6 PUFAs from other FA groups)
## 4.     FA_2: FA grouping 2 (separates out PUFAs from other FA groups)
## 5.     Resource: Resources typically associated with FA
## 6.     Specific_PUFA: Short name for selected PUFAs
## 7.     Long_name: long name of FA

* Missing data codes: NA


Third folder: "03_Dispersal_traits" 
---------------------------------------

* NOTES: 
	* Data and R script for correlation heatmap of selected aquatic invertebrate dispersal traits

Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 3
* Total file size: 36 KB
* Range of individual file sizes: 3 KB - 23 KB
* File formats: .csv, .R, and .txt

Table of Contents
-----------------

"RScript_correlation_invertebrate_dispersal_trait.R"
"Burdon_et_al_Fig.6_invert_dispersal_trait_data.csv"
"README_Burdon_et_al_Invert_dispersal_traits.txt"


Details for:"Burdon_et_al_Fig.6_invert_dispersal_trait_data.csv"
---------------------------------------

* Description: A comma-delimited file containing community-weighted mean abundances for selected aquatic invertebrate dispersal trait modalities.

* Format(s): .csv

* Size(s): 23 KB

* Dimensions: 101 rows by 18 columns

* Variables:
## 1.     Site: Site number	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Site_identifier: Unique code given to each site based on country, site, and site type 
## 4.     Site_code_2: Unique code given to each site in each case-study basin
## 5.     Site_name: Site name	
## 6.     Type_1: Site type (forest, unbuffered, buffered)
## 7.     Type_2: Site type (paired or reference)
## 8.     Site_block: Blocking term for paired sites (unbuffered, buffered)
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

* Missing data codes: NA

Fourth folder: "04_Regression_analyses" 
---------------------------------------

* NOTES: 
	* Data and R scripts for linear mixed effects (LME) models testing correlations between realized trophic connectivity in riparian spiders and their aquatic prey

Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 3
* Total file size: 273 KB
* Range of individual file sizes: 4 KB - 122 KB
* File formats: .csv, .R, and .txt

Table of Contents
-----------------

"RScript_mixed_effects_regression_analyses.R"
"Burdon_et_al_Fig.7_multivariate_data.csv"
"README_Burdon_et_al_LME_regression_analyses.txt"


Details for:"Burdon_et_al_Fig.7_multivariate_data.csv"
---------------------------------------

* Description: A comma-delimited file containing data for riparian spiders, their aquatic prey, and environment properties from stream sites in Sweden, Norway, Romania, and Belgium

* Format(s): .csv

* Size(s): 122 KB

* Dimensions: 251 rows by 38 columns

* Variables:
## 1.     Site_no: Site number	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Taxa_code: Short name for arachnid family or taxonomic group	
## 4.     Genus_family: Family or taxonomic groups used for arachnids	
## 5.     Site_name: Site name
## 6.     Type_1: Site type (forest, unbuffered, buffered)
## 7.     Type_2: Site type (paired or reference)
## 8.     Site_block: Blocking term for paired sites (unbuffered, buffered)
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

* Missing data codes: NA



Fifth folder: "05_SEM"
---------------------------------------

* NOTES:
	* This R script performs the structural equation modelling (SEM). 
	* The main goal is to assess the correlation between realized trophic connectivity in riparian spiders (indicated by fatty acids).
	* The SEMs do this by accounting for key contingencies related to spider community composition, prey aquatic dispersal traits, algal productivity, and environmental context 
	*  Environmental context includes ecosystem size indicated by stream width, riparian condition PC1, and catchment impact PC1
	*  This script also generates the results used in Figures 8 and S8 reported in Burdon et al. (in press) Ecological Monographs
	*  This script also generates the Table 5 results reported in Burdon et al. (in press) Ecological Monographs

Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 5
* Total file size: 476 KB
* Range of individual file sizes: 5 KB - 296 KB
* File formats: .csv, .R, and .txt

Table of Contents
-----------------

"RScript_ecological_connectivity_SEM_PUFAs.R"
"Burdon_et_al_Fig.8_multivariate_data.csv"
"Burdon_et_al_Fig.8_Tetragnathidae_SEM_data_std.csv"
"Burdon_et_al_Fig.8_Lycosidae_SEM_data_std.csv"
"README_Burdon_et_al_SEM_analyses.txt"


Details for:"Burdon_et_al_Fig.8_multivariate_data.csv"
---------------------------------------

* Description: A comma-delimited file containing multivariate data for SEMs including fatty acid composition of riparian spiders, aquatic prey, and environmental variables

* Format(s): .csv

* Size(s): 296 KB

* Dimensions: 257 rows by 69 columns

* Variables:
## 1.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 2.     Taxa_code: Short name for arachnid family or taxonomic group	
## 3.     Type_1: Site type (forest, unbuffered, buffered)
## 4.     Type_2: Site type (paired or reference)
## 5.     Site_block: Blocking term for paired sites (unbuffered, buffered)
## 6.     Site_no: Site number at the country level
## 7.     Site_name: Site name
## 8.     Group: Arachnida	
## 9.     Genus_family: Family or taxonomic groups used for arachnids	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest level: Lowest taxonomic level used for arachnids
## 12.    Tribe
## 13-49. FA1-FA37: Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
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

* Missing data codes: NA	


Details for:"Burdon_et_al_Fig.8_Tetragnathidae_SEM_data_std"
---------------------------------------

* Description: A comma-delimited file containing transformed and standardized multivariate data for SEMs including fatty acid composition of tetragnathid riparian spiders, aquatic prey, and environmental variables

* Format(s): .csv

* Size(s): 63 KB

* Dimensions: 57 rows by 69 columns

* Variables:
## 1.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 2.     Taxa_code: Short name for arachnid family or taxonomic group	
## 3.     Type_1: Site type (forest, unbuffered, buffered)
## 4.     Type_2: Site type (paired or reference)
## 5.     Site_block: Blocking term for paired sites (unbuffered, buffered)
## 6.     Site_no: Site number at the country level
## 7.     Site_name: Site name
## 8.     Group: Arachnida	
## 9.     Genus_family: Family or taxonomic groups used for arachnids	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest level: Lowest taxonomic level used for arachnids
## 12-48. FA1-FA37: Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
## 49.    EPA_perc: % concentration of Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 50.    DHA_perc: % concentration of Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 51.    LNA_perc: % concentration of Linoleic acid (LNA), an omega-6 fatty acid
## 52.    ALA_perc: % concentration of Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 53.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 54.    chl_a_m2_day: Algal accrual on tiles
## 55.    Width: Stream wetted channel width (m)
## 56.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 57.    Catchment_PC1: Catchment impact PC1 (42.3%)	
## 58.     s3: Maximum body size (cm) ≥ 0.5–1
## 59.    cd1: Life cycle duration ≤ 1 year
## 60.    cy2: Potential number of reproductive cycles per year - One
## 61.    dis3: Dispersal strategy Aerial passive
## 62.    dis4: Dispersal strategy Aerial active
## 63.    life1: Adult life span < 1 week
## 64.    life2: Adult life span ≥ 1 week – 1 month
## 65.    wnb2: Wing pair type - 1 pair + halters
## 66.    wnb3: Wing pair type - 1 pair + 1 pair of small hind wings
## 67.    wnb5: Wing pair type - 2 similar-sized pairs
## 68.    EPA_ALA_ratio: Ratio of EPA to ALA
## 69. Tetragnathidae: Tetragnathid spider catch per unit effort (individuals per m2 per min) 

* Missing data codes: NA	

Details for:"Burdon_et_al_Fig.8_Lycosidae_SEM_data_std"
---------------------------------------

* Description: A comma-delimited file containing transformed and standardized multivariate data for SEMs including fatty acid composition of lycosid riparian spiders, aquatic prey, and environmental variables

* Format(s): .csv

* Size(s): 77 KB

* Dimensions: 70 rows by 69 columns

* Variables:
## 1.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 2.     Taxa_code: Short name for arachnid family or taxonomic group	
## 3.     Type_1: Site type (forest, unbuffered, buffered)
## 4.     Type_2: Site type (paired or reference)
## 5.     Site_block: Blocking term for paired sites (unbuffered, buffered)
## 6.     Site_no: Site number at the country level
## 7.     Site_name: Site name
## 8.     Group: Arachnida	
## 9.     Genus_family: Family or taxonomic groups used for arachnids	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest level: Lowest taxonomic level used for arachnids
## 12-48. FA1-FA37: Relative concentrations of 37 fatty acids (FAs) - NOTE: FA38 excluded as not detected in all countries
## 49.    EPA_perc: % concentration of Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 50.    DHA_perc: % concentration of Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 51.    LNA_perc: % concentration of Linoleic acid (LNA), an omega-6 fatty acid
## 52.    ALA_perc: % concentration of Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 53.    Spider_body_size: Community-weighted mean (CWM) abundances for spider body size at sites spiders collected from
## 54.    chl_a_m2_day: Algal accrual on tiles
## 55.    Width: Stream wetted channel width (m)
## 56.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 57.    Catchment_PC1: Catchment impact PC1 (42.3%)	
## 58.     s3: Maximum body size (cm) ≥ 0.5–1
## 59.    cd1: Life cycle duration ≤ 1 year
## 60.    cy2: Potential number of reproductive cycles per year - One
## 61.    dis3: Dispersal strategy Aerial passive
## 62.    dis4: Dispersal strategy Aerial active
## 63.    life1: Adult life span < 1 week
## 64.    life2: Adult life span ≥ 1 week – 1 month
## 65.    wnb2: Wing pair type - 1 pair + halters
## 66.    wnb3: Wing pair type - 1 pair + 1 pair of small hind wings
## 67.    wnb5: Wing pair type - 2 similar-sized pairs
## 68.    EPA_ALA_ratio: Ratio of EPA to ALA
## 69. Lycosidae: Lycosid spider catch per unit effort (individuals per m2 per min) 

* Missing data codes: NA	




Sixth folder: "06_Algae_spider_FAs"
---------------------------------------

* NOTES:
	* This R script performs the generalised additive modelling (GAMs) used to explore non-linear relationship between riparian spiders and stream periphton
	* Specifically I wanted to test if there was a "subsidy-stress" relationships between algal production and realized ecological connectivity


Summary Metrics
---------------

* Major sub-folder count: 0
* File count: 3
* Total file size: 100 KB
* Range of individual file sizes: 3 KB - 50 KB
* File formats: .csv, .R, and .txt

Table of Contents
-----------------

"RScript_algae_spider_FAs_GAMs.R"
"Burdon_et_al_Fig.S7_multivariate_data.csv"
"README_Burdon_et_al_Algae_spider_FAs.txt"
"Burdon_et_al_Fig.S7.png"


Details for:"Burdon_et_al_Fig.8_multivariate_data.csv"
---------------------------------------

* Description: A comma-delimited file containing multivariate data for SEMs including fatty acid composition of riparian spiders, aquatic prey, and environmental variables

* Format(s): .csv

* Size(s): 50 KB

* Dimensions: 257 rows by 22 columns

* Variables:
## 1.     Genus_family: Family or taxonomic groups used for arachnids	
## 2.     Country: Country where spider was collected (Sweden, Norway, Romania, Belgium)	
## 3.     Taxa_code: Short name for arachnid family or taxonomic group	
## 4.     Type_1: Site type (forest, unbuffered, buffered)
## 5.     Type_2: Site type (paired or reference)
## 6.     Site_block:	Blocking term for paired sites (unbuffered, buffered)
## 7.     Site_no: Site number	
## 8.     Site_name: Site name	
## 9.     Group: Arachnida	
## 10.    Mode: Mode of hunting for arachnids
## 11.    Lowest_level: Lowest taxonomic level used	
## 12.    Tribe: Used for ground beetles not included in this dataset
## 13.    Ratio_DHA_LNA: Ratio of DHA to LNA
## 14.    Ratio_EPA_LNA: Ratio of EPA to LNA
## 15.    Ratio_EPA_ALA: Ratio of EPA to ALA
## 16.    EPA_perc: %EPA, Eicosapentaenoic acid (EPA), an omega-3 fatty acid 
## 17.    DHA_perc: %DHA, Docosahexaenoic acid (DHA), an omega-3 fatty acid 
## 18.    LNA_perc: %LNA, Linoleic acid (LNA), an omega-6 fatty acid
## 19.    ALA_perc: #ALA, Alpha-linolenic acid (ALA), an omega-3 fatty acid 
## 20.    chl_a_m2_day: Algal accrual on tiles
## 21.    Riparian_PC1: Riparian condition PC1 (36.7%)
## 22     Catchment_PC1: Catchment impact PC1 (42.3%)	pairs
## 69.    EPA_ALA_ratio: Ratio of EPA to ALA

* Missing data codes: NA	

