## Factors influencing the biomass of large-bodied parrotfishes in the absence of fishing on coral reefs in Florida, USA

This repository includes data and analysis scripts for the following article:

Zuercher, R., D. Kochan, A. Harborne. Factors influencing the biomass of large-bodies parrotfishes in the absence of fishing on coral reefs in Florida, USA. *Submitted to Coral Reefs*

---
## Collaborators:
Alastair Harborne, Tropical Fish Ecology Lab at Florida International University (FIU)
Rachel Zuercher, Tropical Fish Ecology Lab at FIU + University of Rhode Island (URI) Coastal Resources Center
David Kochan, Tropical Fish Ecology Lab at FIU + Florida Fish and Wildlife Conservation Commission   

**Contacts**: rachel.zuercher@gmail.com; aharborne@fiu.edu

---
## Description:
Parrotfishes are a functionally critical component of Caribbean reef fish assemblages, with large-bodied parrotfish species exerting particularly important top-down control on macroalgae. Despite their importance, low biomasses of large-bodied parrotfishes on many reefs hampers our ability to study and understand their ecology. Florida reefs, where most parrotfish fishing has been illegal since 1992, present a unique opportunity to quantify the variables influencing their biomass. Using boosted regression tree models, we identify the major predictors of four species of Atlantic large-bodied parrotfish. Positive relationships between parrotfish presence and reef complexity, the area of surrounding reef, and availability of seagrass habitat identify these as key variables for the ecology of large-bodied parrotfishes. Strong positive relationships between parrotfish presence and biomass and the biomass of other parrotfishes on a reef suggest little intraspecies competition and that all four species are responding to a similar set of drivers. However, relationships between parrotfish presence and biomass and depth, coral cover and the proximity of a reef to deep water habitats differed markedly among species, highlighting distinct habitat preferences among the large-bodied parrotfishes. The biomass of Stoplight Parrotfish showed a significant negative correlation with human population density, alluding to negative impacts of anthropogenic stressors on this species. These results improve our ability to target important biophysical drivers of large-bodied parrotfishes with management approaches

--- 
## In the repository:
Four .csv files and one .R script are needed to replicate these analyses.

`RVC_impact.csv` -- data file containing each Reef Visual Census survey site used in the fishing impact model (rows), the snapper-grouper biomass for each site (column), and all explanatory variables considered for the fishing impact model (columns)

`ReefPoints_impact.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the fishing impact model

`RVC_biomass.csv` -- data file containing each Reef Visual Census survey site used in the fish biomass models (rows), biomass of all species groups (columns), and all explanatory variables considered for the biomass models, including fishing impact estimated by this project and extrapolated to these sites in ArcGIS (columns)

`ReefPoints_biomass.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the biomass models

`MOW_FL_3.2021.R` script -- runs all analyses and creates plots for Zuercher et al. 2021

The analyses can be replicated by changing the working directory in `MOW_FL_3.2021.R` to the location on your computer where you have stored the .R and .csv files. Additional analyses for this project were conducted in ArcGIS Pro. Spatial data layers are housed privately, but can be requested for the purpose of replication or for additional research. Questions about the code and requests for spatial data layers should be directed to Rachel Zuercher (rachel.zuercher@gmail.com).
