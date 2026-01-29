library(BayesMaxEnt)
library(maxnet)
library(disdat)
# library(bayesplot)
library(parallel)

source("MaxEntBayesFunctions.R")

RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")



FitAll(dataname="AWT", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="CAN", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="NSW", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="NZ", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="SA", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="SWI", remove=RemoveNames, classes="l", validate = FALSE)

FitAll(dataname="AWT", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="CAN", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="NSW", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="NZ", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="SA", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="SWI", remove=RemoveNames, classes="l", fit = FALSE)


# "AWT", 
#. Needs to specify group
EnvPO <- names(disBg("AWT"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- unique(disPo("AWT")$spid)

# FitAndValidate(sp=AWTSpecies[1], DataName="AWT", Env=EnvPO, nclust=1, small=TRUE)

# AWTThing <- sapply(AWTSpecies[1:3], FitAndValidate, DataName="AWT", Env=EnvPO, nclust=6, small=TRUE)
# AWTThing <- sapply(AWTSpecies, FitAndValidate, DataName="AWT", Env=EnvPO, nclust=6)
AWTThing <- sapply(Species, JustFit, DataName="AWT", Env=EnvPO, classes="l")
AWTThing <- sapply(Species, JustValidate, DataName="AWT", Env=EnvPO, small=FALSE, 
                  nclust = 6, classes="l", verbose=TRUE)


# "CAN"
EnvPO <- names(disBg("CAN"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- unique(disPo("CAN")$spid)
# CanThing <- sapply(CanSpecies, FitAndValidate, DataName="CAN", Env=EnvPO, nclust=6)

CanThing <- sapply(Species, JustFit, DataName="CAN", Env=EnvPO, classes="l")
CanThing <- sapply(Species, JustValidate, DataName="CAN", Env=EnvPO, small=FALSE, 
                   nclust = 6, classes="l", verbose=TRUE)


# "NSW", 
#. Needs to specify group
EnvPO <- names(disBg("NSW"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- unique(disPo("NSW")$spid)
# NSWThing <- sapply(NSWSpecies, FitAndValidate, DataName="NSW", Env=EnvPO, nclust=6, small=TRUE)

NSWThing <- sapply(Species, JustFit, DataName="NSW", Env=EnvPO, classes="l")
NSWThing <- sapply(Species, JustValidate, DataName="NSW", Env=EnvPO, small=FALSE, 
                   nclust = 6, classes="l", verbose=TRUE)

# "NZ", 
EnvPO <- names(disBg("NZ"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- unique(disPo("NZ")$spid)

# NZThing <- sapply(NZSpecies, JustFit, DataName="NZ", Env=EnvPO, classes="l", nclust=6)

# JustFit(sp=NZSpecies[1], DataName="NZ", Env=EnvPO, small=TRUE, classes="l", verbose=TRUE)
# JustValidate(sp=NZSpecies[1], DataName="NZ", Env=EnvPO, small=FALSE, nclust = 6, classes="l", verbose=TRUE)

NZThing <- sapply(Species, JustFit, DataName="NZ", Env=EnvPO, classes="l")
NZThing <- sapply(Species, JustValidate, DataName="NZ", Env=EnvPO, small=FALSE, 
                  nclust = 6, classes="l", verbose=TRUE)

# thing <- FitAndValidate(sp=NZSpecies[1], DataName="NZ", Env=EnvPO, nclust=6, classes="lq", small = TRUE)


# "SA", 
EnvPO <- names(disBg("SA"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
# Species <- "can02"
Species <- unique(disPo("SA")$spid)

SAThing <- sapply(SASpecies, FitAndValidate, DataName="SA", Env=EnvPO, nclust=6)

AWTThing <- sapply(AWTSpecies, JustFit, DataName="NZ", Env=EnvPO, classes="l")
AWTThing <- sapply(AWTSpecies, JustValidate, DataName="NZ", Env=EnvPO, small=FALSE, 
                   nclust = 6, classes="l", verbose=TRUE)


# "SWI"
EnvPO <- names(disBg("SWI"))
EnvPO <- EnvPO[!(EnvPO%in%RemoveNames)]
Species <- "can02"
SWISpecies <- unique(disPo("SWI")$spid)

SWIthing <- sapply(SWISpecies, FitAndValidate, DataName="SWI", Env=EnvPO, nclust=6)

AWTThing <- sapply(AWTSpecies, JustFit, DataName="NZ", Env=EnvPO, classes="l")
AWTThing <- sapply(AWTSpecies, JustValidate, DataName="NZ", Env=EnvPO, small=FALSE, 
                   nclust = 6, classes="l", verbose=TRUE)


# mcmc_trace(BigFit$mcmc, pars = c("alpha", "beta[1]"))
