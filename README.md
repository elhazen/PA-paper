# PA-paper

Code supporting the publication "Where do they not go?". 

Code authors: Elliott Hazen (NOAA, UCSC), Heather Welch (UCSC, NOAA), Stephanie Brodie (UCSC, NOAA)

Relevant manuscripts:

Hazen et al. 2018 “A dynamic ocean management tool to reduce bycatch and support sustainable fisheries.” Science Advances 4: eaar3001.

Welch et al. 2018 "Practical considerations for operationalizing dynamic management tools." Journal of Applied Ecology.DOI: 10.1111/1365-2664.13281.

Abrahms et al. 2019 "WhaleWatch: a dynamic management tool for predicting blue whale density in the California Current." Journal of Applied Ecology 54: 1415-1428.

Tsalyuk et al. 2019 "" Ecological Monographs

Description of scripts:
RunAllCode.R - Master file to run all the relevant code assuming you have tracking data.
load_packages_PA.R - Loads all relevant packages.

2 - generate-pseudo-absences.R Creates 4 pseudoabsence types: correlated random walks, reverse correlated random walks, buffer absences, background absences for both blue whales and elephants

3 - Extract-BW.R & 3 - Extract-Elephant.R

Samples downloaded environmental data as rasterstacks - downloaded from the data provider, such as SSHa via AVISO / CMEMS, e.g. https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/sla-h.html. Elephant data were static and sampled as shapefiles

4 - BuildModels.R - Fits Generalized Linear Mixed Models, Generalized Additive Mixed Models, and Boosted Regression Tree models for both Blue Whale and Elephant species.

5 - EvaluateModels.R - Calculates evaluation statistics for 100%, leave-one-out by space, leave-one-out by time (for blue whales)

6 - FinalFigs.R - Code to support creating final figures for the manuscript 
