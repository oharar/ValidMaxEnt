library(BayesMaxEnt)
library(maxnet)
library(disdat)
# library(bayesplot)
library(parallel)

source("MaxEntBayesFunctions.R")

RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")

# First, fit MaxEnt models to PO data
# Code for each data set, because in practice we don't want to run them all together
# Note that we could run these by species, using WriteBashScript() to write a script, 
# but these are quick enough to run together

FitAll(dataname="AWT", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="CAN", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="NSW", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="NZ", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="SA", remove=RemoveNames, classes="l", validate = FALSE)
FitAll(dataname="SWI", remove=RemoveNames, classes="l", validate = FALSE)



# Now validate the models by using the output to calculate the linear predictors,
#  and use these in a logistic regression to model the PO data
FitAll(dataname="AWT", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="CAN", remove=RemoveNames, classes="l", fit = FALSE, overwrite = TRUE)
FitAll(dataname="NSW", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="NZ", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="SA", remove=RemoveNames, classes="l", fit = FALSE)
FitAll(dataname="SWI", remove=RemoveNames, classes="l", fit = FALSE)


