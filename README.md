# PA-paper

Code supporting the publication "Where do they not go?". 

Tracking data offer where species choose to go but not where they don't. Generating pseudo-absences is an approach to simulating where they could have gone but did not go, yet the impacts of such generation have not yet been tested. The goal of this repo is to assimilate a suite of environmental data along with relevant tracking data for generating absences, sampling environmental data, model fitting, and prediction and evaluation for the manuscript. Ultimately, we find that traditional model predictive skill metrics may not identify the best model depending on the goals and desires when building a predictive model, and we reiterate the need to focus on biological realism when evaluating such models.

Code authors: Elliott Hazen (NOAA, UCSC), Heather Welch (UCSC, NOAA), Stephanie Brodie (UCSC, NOAA), Briana Abrahms (UW), Gemma Carroll (UCSC, NOAA)

Relevant manuscripts:

Abrahms, B., Welch, H., Brodie, S., Jacox, M.G., Becker, E.A., Bograd, S.J., Irvine, L.M., Palacios, D.M., Mate, B.R. and Hazen, E.L., 2019. Dynamic ensemble models to predict distributions and anthropogenic risk exposure for highly mobile species. Diversity and Distributions, 25(8), pp.1182-1193.

Tsalyuk, M., Kilian, W., Reineking, B. and Getz, W.M., 2019. Temporal variation in resource selection of African elephants follows long‐term variability in resource availability. Ecological Monographs, 89(2), p.e01348.

Hazen, E.L., Scales, K.L., Maxwell, S.M., Briscoe, D.K., Welch, H., Bograd, S.J., Bailey, H., Benson, S.R., Eguchi, T., Dewar, H. and Kohin, S., 2018. A dynamic ocean management tool to reduce bycatch and support sustainable fisheries. Science advances, 4(5), p.eaar3001.

Description of scripts:

RunAllCode.R - Master file to run all the relevant code assuming you have tracking data.

load_packages_PA.R - Loads all relevant packages.
PseudoFunctions.R - Functions to create 4 types of pseudotracks.
get_EnvData.R - Functions to assimilate relevant environmental data for model fitting and evaluation
ExtractionFunction.R - Functions to sample downloaded environmental data.
ModelFunctions.R - Functions to build GLMM, GAMM, BRT.
Model_Eval_Fcns.R - Functions to evaluate GLMM, GAMM, BRT models.


2 - generate-pseudo-absences.R Creates 4 pseudoabsence types: correlated random walks, reverse correlated random walks, buffer absences, background absences for both blue whales and elephants

3 - Extract-BW.R & 3 - Extract-Elephant.R

Samples downloaded environmental data as rasterstacks - downloaded from the data provider, such as SSHa via AVISO / CMEMS, e.g. https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/sla-h.html. Elephant data were static and sampled as shapefiles

4 - BuildModels.R - Fits Generalized Linear Mixed Models, Generalized Additive Mixed Models, and Boosted Regression Tree models for both Blue Whale and Elephant species.

5 - EvaluateModels.R - Calculates evaluation statistics for 100%, leave-one-out by space, leave-one-out by time (for blue whales)

6 - FinalFigs.R - Code to support creating final figures for the manuscript 
