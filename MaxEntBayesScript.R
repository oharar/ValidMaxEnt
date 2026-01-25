# First argument is data set, no others used yet.
# second is classes: subset of "lqpht" 
args <- commandArgs(trailingOnly = FALSE)

dataname <- args[1]
if(length(args)>1) {
 classes <- args[2] 
} else (
  classes <- "l"
)
library(BayesMaxEnt)
library(maxnet)
library(disdat)
library(bayesplot)
library(parallel)

source("MaxEntBayesFunctions.R")
RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")
dataname <- args[1]

# Fit...
Env <- names(disBg(dataname))
Species <- unique(disPo(dataname)$spid)
Thing <- sapply(Species, FitAndValidate, DataName=dataname, Env=Env, nclust=6, small=TRUE)

