## Factors influencing the biomass of large-bodied parrotfishes in the absence of fishing on coral reefs in Florida, USA

This repository includes data and analysis scripts for the following article:

Zuercher, R., Kochan, D., & Harborne,A. R. (2023). Factors influencing the biomass of large-bodiedparrotfish species in the absence of fishing on coral reefs inFlorida, USA.Journal of Fish Biology,1–12.https://doi.org/10.1111/jfb.15557

---
## Collaborators:
*Alastair Harborne, Florida International University*

*Rachel Zuercher, Florida International University; University of Rhode Island Coastal Resources Center*

*David Kochan, Florida International University; Florida Fish and Wildlife Conservation Commission*   

**Contacts**: rachel.zuercher@gmail.com; aharborne@fiu.edu

---
## Description:
Parrotfishes are a functionally critical component of Caribbean reef fish assemblages,with large-bodied parrotfish species exerting particularly important top-down controlon macroalgae. Despite their importance, low biomasses of large-bodied parrotfisheson many reefs hamper our ability to study and understand their ecology. Floridareefs, where most parrotfish fishing has been illegal since 1992, present a uniqueopportunity to explore covariates of their distribution. Using boosted regression treemodels and 23 covariates, this study identified the major predictors of four speciesof Atlantic large-bodied parrotfishes. Maximum hard substrate relief, the area of thesurrounding reef, and the availability of seagrass habitat were each positively relatedto parrotfish presence. Strong positive relationships between parrotfish presence andbiomass and the biomass of other parrotfishes on a reef suggest that all four speciesresponded to a similar subset of environmental conditions. However, relationshipsbetween parrotfish presence and biomass and depth, habitat type, coral cover, andthe proximity of a reef to deepwater habitats differed among species, highlightingdistinct habitat preferences. These results can improve managers’ability to targetimportant biophysical correlates of large-bodied parrotfishes with appropriate man-agement interventions and identify areas for protection.

--- 
## In the repository:
One .csv files and one .R script are needed to replicate these analyses.

`RVC_parrotfish_2012-2018.csv` -- data file containing each Reef Visual Census survey site used in the parrotfish presence/absence and biomass models, a column for the presence or absence (0 or 1) of each of four species: Midnight, Blue, Rainbow and Stoplight Parrotfishes, a column for the biomass of each of the four focal species (g/m2), and a column for each explanatory variable considered in the models


`MOW_parrotfish_4.2022.R` script -- runs all analyses and creates plots for Zuercher et al. "Factors influencing the biomass of large-bodied parrotfishes..."

The analyses can be replicated by changing the working directory in `MOW_parrotfish_4.2022.R` to the location on your computer where you have stored the .R and .csv files. Spatial data layers used for this research are housed privately, but can be requested for the purpose of replication or for additional research. Questions about the code and requests for spatial data layers should be directed to Rachel Zuercher (rachel.zuercher@gmail.com).
