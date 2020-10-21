### running model eval functions
### written by HW, SB, BA 2019

### load packages ####

soure("Model_Eval_Fcns.R")

### define global terms ####
#see description of terms in Model_Eval_Fcns.R
gbm.x.bwhale <- c("Bathymetry","Oxygen_100m", "SST", "Chla_25km_monthly", 
           "SLA", "Oxygen_Surface", "SST_SD", "Chla_4km_8day", 
           "MLD", "Rugosity","U", "FSLE_max", "eke", "V",
           "Theta_max")

gbm.x.ephant <- c()

lr=0.005
tc=5
bf=.75
nt = 2000
