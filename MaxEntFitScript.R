# Usage
# Rscript MaxEntFitScript.R --args NZ 
# Rscript MaxEntFitScript.R --args NZ l

# First argument is data set (One of "AWT", "CAN", "NSW", "NZ", "SA", "SWI")
# second is classes: subset of "lqpht" 

args <- commandArgs(trailingOnly = TRUE)
dataname <- as.character(args[2])

# cat(dataname, "\n", file="Try.txt")

if(length(args)>2) {
 classes <- args[3] 
} else (
  classes <- "l"
)
 cat(dataname, "\n", classes, "\n", file="Try.txt")

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


 # This is for testing
#  Thing <- sapply(Species[1:3], JustFit, DataName=dataname, Env=Env, small=TRUE, classes=classes)
# This is the real thing 
 Thing <- sapply(Species, JustFit, DataName=dataname, Env=Env, classes=classes)
