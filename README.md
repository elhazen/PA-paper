# PA-paper

Code supporting the publication "Where do they not go?". 

Tracking data offer where species choose to go but not where they don't. Generating pseudo-absences is an approach to simulating where they could have gone but did not go, yet the impacts of such generation on species distribution model fit and performance have not yet been tested. The goal of this repo is to use tracking data to generate a suite of pseudo-absences, sample environmental data (using , and to fit, predict, and evaluate species distribution models. Ultimately, we find that traditional model performance metrics alone may not identify the best model depending on the goals and desires when building a predictive model, and we reiterate the need to focus on biological realism along-side traditional skill metrics when evaluating such models.

EDIT 03/31/2023: Consider using ANIMOTUM package for pseudotrack simulation from SSM parameters rather than just speed and turn angle (this repo). Code here: https://ianjonsen.github.io/aniMotum/index.html and paper here: https://doi.org/10.1111/2041-210X.14060

**Code authors:** Elliott Hazen (NOAA, UCSC), Heather Welch (UCSC, NOAA), Stephanie Brodie (UCSC, NOAA), Briana Abrahms (UW), Gemma Carroll (UCSC, NOAA)

**Relevant manuscripts:**

Abrahms, B., Welch, H., Brodie, S., Jacox, M.G., Becker, E.A., Bograd, S.J., Irvine, L.M., Palacios, D.M., Mate, B.R. and Hazen, E.L., 2019. Dynamic ensemble models to predict distributions and anthropogenic risk exposure for highly mobile species. Diversity and Distributions, 25(8), pp.1182-1193.

Tsalyuk, M., Kilian, W., Reineking, B. and Getz, W.M., 2019. Temporal variation in resource selection of African elephants follows long‐term variability in resource availability. Ecological Monographs, 89(2), p.e01348.

Hazen, E.L., Scales, K.L., Maxwell, S.M., Briscoe, D.K., Welch, H., Bograd, S.J., Bailey, H., Benson, S.R., Eguchi, T., Dewar, H. and Kohin, S., 2018. A dynamic ocean management tool to reduce bycatch and support sustainable fisheries. Science advances, 4(5), p.eaar3001.

**Required data: **

I) Presence data (tracks) with lat, lon, time, and ideally an error estimate for positions to inform environmental extraction. Data are expected to have already been quality controlled for erroneous points. 
II) Environmental data for blue whales (rasters) can be downloaded from multiple sources. The NetCDF files have been converted to a regularized grid for sampling by presences and pseudoabsences before extraction in this example. Other approaches to environmental download and extraction can be used (e.g. 2, 3a, 3b - https://github.com/elhazen/EcoCast-SciAdv).
For elephants, environmental layers (shapefiles) are used to calculate distance to roads, distance to water, and long-term average NDVI. 

**Description of scripts:**

RunAllCode.R - Master file to run all the relevant code assuming you have tracking data.

**Functions:** 

load_packages_PA.R - Loads all relevant packages.
PseudoFunctions.R - Functions to create 4 types of pseudoabsences (described below).
get_EnvData.R - Functions to assimilate relevant environmental data for model fitting and evaluation.
ExtractionFunction.R - Functions to sample downloaded environmental data at tagging and pseudo-absence locations.
ModelFunctions.R - Functions to build GLMM, GAMM, BRT.
Model_Eval_Fcns.R - Functions to evaluate GLMM, GAMM, BRT models.

**Scripts that call functions:** 

2 - generate-pseudo-absences.R Creates four pseudoabsence types: correlated random walks, reverse correlated random walks, buffer absences, background absences for both blue whales and elephants. Calls PseudoFunctions.R.

3 - Extract-BW.R Extracts environmental data to tagging data and pseudo-absences. Samples downloaded environmental data as rasterstacks - downloaded from the data provider, such as SSHa via AVISO / CMEMS, e.g. https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/sla-h.html. 

3 - Extract-Elephant.R - Elephant data were static and sampled using shapefiles.

4 - BuildModels.R - Fits Generalized Linear Mixed Models, Generalized Additive Mixed Models, and Boosted Regression Tree models for both Blue Whale and Elephant species. Calls ModelFunctions.R

5 - EvaluateModels.R - Calculates model performanceevaluation metrics such as AUC, TSS, explained deviance via multiple cross-validation methods:statistics for 100%, leave-one-out by space, leave-one-out by time (for blue whales). Calls Model_Eval_Fcns.R

6 - FinalFigs.R - Code to support creating final figures for the manuscript. 

