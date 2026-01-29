# Usage
# RScript MaxEntBayesScript.R --args NZ 
# RScript MaxEntBayesScriptSpecies.R --args NZ nz01 l

lqpht

# First argument is data set (One of "AWT", "CAN", "NSW", "NZ", "SA", "SWI")
# second argument is species
# third is classes: subset of "lqpht" 

args <- commandArgs(trailingOnly = TRUE)
dataname <- as.character(args[2])
species <- as.character(args[3])

# cat(dataname, "\n", file="Try.txt")

if(length(args)>3) {
 classes <- args[4] 
} else (
  classes <- "l"
)
# cat(dataname, "\n", classes, "\n", file="Try.txt")

library(BayesMaxEnt)
library(maxnet)
library(disdat)
library(parallel)

source("MaxEntBayesFunctions.R")
RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")

# Fit...
Env <- names(disBg(dataname))
Env <- Env[!(Env%in%RemoveNames)]

# Species <- unique(disPo(dataname)$spid)
# Thing <- sapply(Species, FitAndValidate, DataName=dataname, Env=Env, nclust=6, classes=classes)

 # This is for testing
#Thing <- sapply(Species[1:3], JustFit, DataName=dataname, Env=Env, small=TRUE)
 Thing <- JustFit(species, DataName=dataname, Env=Env, small=TRUE)
