## Factors influencing the biomass of large-bodied parrotfishes in the absence of fishing on coral reefs in Florida, USA

This repository includes data and analysis scripts for the following article:

Zuercher, R., D. Kochan, A. Harborne. Factors influencing the biomass of large-bodies parrotfishes in the absence of fishing on coral reefs in Florida, USA. *Submitted to Coral Reefs*

---
## Collaborators:
*Alastair Harborne, Florida International University*

*Rachel Zuercher, Florida International University; University of Rhode Island Coastal Resources Center*

*David Kochan, Florida International University; Florida Fish and Wildlife Conservation Commission*   

**Contacts**: rachel.zuercher@gmail.com; aharborne@fiu.edu

---
## Description:
Parrotfishes are a functionally critical component of Caribbean reef fish assemblages, with large-bodied parrotfish species exerting particularly important top-down control on macroalgae. Despite their importance, low biomasses of large-bodied parrotfishes on many reefs hampers our ability to study and understand their ecology. Florida reefs, where most parrotfish fishing has been illegal since 1992, present a unique opportunity to quantify the variables influencing their distribution. Using boosted regression tree models and 23 covariates, we identify the major predictors of four species of Atlantic large-bodied parrotfishes. Positive relationships between parrotfish presence and reef complexity, the area of surrounding reef, and availability of seagrass habitat identify these as key variables for the ecology of large-bodied parrotfishes. Strong positive relationships between parrotfish presence and biomass and the biomass of other parrotfishes on a reef suggest little interspecies competition and that all four species are responding to a similar subset of drivers. However, relationships between parrotfish presence and biomass and depth, habitat type, coral cover and the proximity of a reef to deep water habitats differed among species, highlighting some distinct habitat preferences among the large-bodied parrotfishes. These results improve our ability to target important biophysical drivers of large-bodied parrotfishes with appropriate management interventions.

--- 
## In the repository:
One .csv files and one .R script are needed to replicate these analyses.

`RVC_parrotfish_2012-2018.csv` -- data file containing each Reef Visual Census survey site used in the parrotfish presence/absence and biomass models, a column for the presence or absence (0 or 1) of each of four species: Midnight, Blue, Rainbow and Stoplight Parrotfishes, a column for the biomass of each of the four focal species (g/m2), and a column for each explanatory variable considered in the models


`MOW_parrotfish_2.2022.R` script -- runs all analyses and creates plots for Zuercher et al. "Factors influencing the biomass of large-bodied parrotfishes..."

The analyses can be replicated by changing the working directory in `MOW_parrotfish_2.2022.R` to the location on your computer where you have stored the .R and .csv files. Spatial data layers used for this research are housed privately, but can be requested for the purpose of replication or for additional research. Questions about the code and requests for spatial data layers should be directed to Rachel Zuercher (rachel.zuercher@gmail.com).
