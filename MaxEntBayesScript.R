# Usage
# RScript MaxEntBayesScript.R --args NZ 
# RScript MaxEntBayesScript.R --args NZ l

# First argument is data set (One of "AWT", "CAN", "NSW", "NZ", "SA", "SWI")
# second is classes: subset of "lqpht" 

args <- commandArgs(trailingOnly = TRUE)

dataname <- as.character(args[2])

cat(dataname, "\n", file="Try.txt")

if(length(args)>2) {
 classes <- args[3] 
} else (
  classes <- "l"
)
library(BayesMaxEnt)
library(maxnet)
library(disdat)
library(parallel)

source("MaxEntBayesFunctions.R")
RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")

# Fit...
Env <- names(disBg(dataname))
Env <- Env[!(Env%in%RemoveNames)]

Species <- unique(disPo(dataname)$spid)
 Thing <- sapply(Species, FitAndValidate, DataName=dataname, Env=Env, nclust=6, classes=classes)

 # This is for testing
 # Thing <- sapply(Species[1:3], FitAndValidate, DataName=dataname, Env=Env, nclust=6, small=TRUE)
 