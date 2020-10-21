# Code library for Pseudoabsence evaluation manuscript for blue whales and elephants. Stages are to:
# 1) source all relevant functions and packages
# 2) generate pseudo-absences for tracking data
# 3) extract environmental data for tracking data and pseudoabsences
# 4) Fit models
# 5) Evaluate models for predictive skill and biological realism
# 6) Make figures for manuscript

# Code written by Elliott Hazen, Heather Welch, Steph Brodie, and Briana Abrahms 
# Date: 10/21/2020

# 1) Source functions and packages
source("LoadPackages_PA.R")
source("PseudoFunctions.R")
source("get_EnvData.R")
source("ExtractionFunction.R")
source("ModelFunctions.R")
source("Model_Eval_Fncs.R")

# 2) generate pseudo-absences for tracking data
### Start with Blue Whales, then run elephants

source("2 - Generate-pseudo-absences.R")

# 3) Extract environmental data

source("3 - Extract-BW.R")
source("3 - Extract-Elephant.R")

# 4) Fit models

source("4 - BuildModels.R")

# 4b) Predict models
source("4b - PredictModels.R")

# 5) Evaluate for model skill and environmental realism

source("5 - EvaluateModels.R")

# 6) Create final figures for manuscript

source("6 - FinalFigs.R")



