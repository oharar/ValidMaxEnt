library(BayesMaxEnt)
library(maxnet)
library(disdat)
library(bayesplot)
library(parallel)

source("MaxEntBayesFunctions.R")

RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")

# "AWT", 
#. Needs to specify group
EnvPO <- names(disBg("AWT"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
AWTSpecies <- unique(disPo("AWT")$spid)

# FitAndValidate(sp=AWTSpecies[1], DataName="AWT", Env=EnvPO, nclust=1, small=TRUE)

AWTThing <- sapply(AWTSpecies[1:3], FitAndValidate, DataName="AWT", Env=EnvPO, nclust=6, small=TRUE)


# "CAN"
EnvPO <- names(disBg("CAN"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
CanSpecies <- unique(disPo("CAN")$spid)
CanThing <- sapply(CanSpecies, FitAndValidate, DataName="CAN", Env=EnvPO, nclust=6)


# "NSW", 
#. Needs to specify group
EnvPO <- names(disBg("NSW"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
NSWSpecies <- unique(disPo("NSW")$spid)
NSWThing <- sapply(NSWSpecies, FitAndValidate, DataName="NSW", Env=EnvPO, nclust=6, small=TRUE)


# "NZ", 
EnvPO <- names(disBg("NZ"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
NZSpecies <- unique(disPo("NZ")$spid)

NZThing <- sapply(NZSpecies, FitAndValidate, DataName="NZ", Env=EnvPO, nclust=6)


# "SA", 
EnvPO <- names(disBg("SA"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- "can02"
SASpecies <- unique(disPo("SA")$spid)

SAThing <- sapply(SASpecies, FitAndValidate, DataName="SA", Env=EnvPO, nclust=6)



# "SWI"
EnvPO <- names(disBg("SWI"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- "can02"
SWISpecies <- unique(disPo("SWI")$spid)

SWIthing <- sapply(SWISpecies, FitAndValidate, DataName="SWI", Env=EnvPO, nclust=6)



# mcmc_trace(BigFit$mcmc, pars = c("alpha", "beta[1]"))
