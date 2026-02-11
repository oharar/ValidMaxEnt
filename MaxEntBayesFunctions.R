
# Plot a species' data
PlotMap <- function(sp, dataname) {
  require(sf)
  require(disdat)
  bg <- disBg(dataname)
  po <- disPo(dataname)
  
  if (dataname %in% c("AWT", "NSW")) {
    grp <- po[po$spid == sp, "group"][1]
    pa <- disPa(dataname, grp)
  } else {
    pa <- disPa(dataname)
  }
  
  plot(sf::st_geometry(disBorder(dataname)), main=sp)
  points(bg$x, bg$y, col="grey80", cex=0.2, pch=16)
  points(pa$x[pa[,sp]==0], pa$y[pa[,sp]==0], col="pink", cex=0.3, pch=16)
  points(pa$x[pa[,sp]==1], pa$y[pa[,sp]==1], col="red", cex=0.5, pch=16)
  points(po$x[po$spid==sp], po$y[po$spid==sp], col="blue", cex=0.5, pch=16)
}



# Function to fit Bayesian MaxEnt model to a sopecies in a DisDat dataset
# species: species name/code
# data: DisDat data
# Envnames: vector of names of environmental variables
# lambda: Value of lambda to use. Defaults to NULL, in which case the value from MaxEnt it used.
# nburn: Number of burn-in iterations, defaults to 10
# niter: Number of iterations after burn-in, defaults to 10
# nchain: Number of chains, defaults to 1

FitModelToDisDat <- function(species, dataname, Envnames, lambda=NULL, 
                             nburn =10, niter = 10, nchain=1, thin=1, classes="l") {
  require(disdat)
  require(maxnet)
  data <- disData(dataname)
  
  SpData <- rbind(data$po[disPo(dataname)$spid==species,c("occ", Envnames)], data$bg[,c("occ", Envnames)])
  EnvScales <- apply(SpData[,Envnames], 2, function(x) c(mean=mean(x), sd=sd(x))) # scale covariates
  SpData[,Envnames] <- apply(SpData[,Envnames], 2, scale) # scale covariates
  
  ToNimble <- SetUpMaxEnt(p=SpData$occ, data=SpData[,Envnames], 
                          regmult = 1, lambda=lambda, classes=classes, 
                          addsamplestobackground=FALSE)
  
  # Fit with MaxNet. This will give us the optimum lambda
  MaxNet.mod <- maxnet(p=ToNimble$Data$y, data=SpData[,Envnames], 
                       f=maxnet.formula(p=ToNimble$Data$y, data=SpData[,Envnames], 
                                        classes=classes))
  if(is.null(lambda)) ToNimble$Const$MeanLambda <- MaxNet.mod$lambda[200]
  
  output <- FitMaxEnt(maxdat=ToNimble, adaptInterval=nburn,
                      nchains=nchain, nburnin = nburn, niter=nburn+niter, thin=thin)
  list(maxnet = MaxNet.mod, mcmc = output)
}



# Fit model
# This fits the MaxEnt (Bayes and classical) to one species
# Arguments:
#  sp: species name
#  DataName: data name, string
#  Env: vector of names of environmental variables
#  small: If TRUE will run 2 short chains, for testing. Defaults to FALSE
#  classes: Classes of MaxEnt features to use. Defaluts to "l", i.e. just linear
#  verbose: Should the function tell us what it has done? Defaults to FALSE
#  overwrite: Should we over-write an old file? Defaults to FALSE.
#  fileprefix: Prefix (incl. folders) for in & output. Defaults to "Results/MaxEntRes" for compatability
JustFit <- function(sp, DataName, Env, small=FALSE, classes="l", verbose=FALSE, overwrite=FALSE,
                    fileprefix="Results/MaxEntRes", lambda=NULL) {
  FileName <- paste0(fileprefix, "MCMC", DataName, "_", classes, "_", sp, ".RData")
  if(small) {
    NBurn <- 10
    NIter <- 10
    NChain <- 2
    Thin=1
    FileName <- gsub("MCMC", "MCMCsmall", FileName)
  } else {
    NBurn <- 1.5e4
    NIter <- 1e3
    NChain <- 5
    Thin <- 5
  }
  if(!file.exists(FileName)|overwrite) {
    if(verbose) message("Started ", sp)
    Fit <- FitModelToDisDat(species=sp, dataname=DataName, Envnames=Env, classes=classes,
                            nburn =NBurn, niter = NIter, nchain=NChain, thin=Thin, lambda=lambda)
    if(verbose) message("Fitted ", sp)
    
    save(Fit, file=FileName)
    if(verbose) message("Done ", sp)
  } else {
    if(verbose) message("Already Done ", sp)
  }
}



# Function to fit calibration model to one MCMC draw
# num.thread=1 used because we are parallelising elsewhere
FitCalib <- function(mcmcpars, env, pres) {
  LP <- env%*%mcmcpars
  LP.c <- LP - mean(LP)
  dat <- data.frame(Pres = pres, LP = LP.c)
  dat$Pred <- dat$LP - mean(dat$LP)
  mod <- INLA::inla(Pres ~ Pred + offset(Pred), data=dat, family = "binomial", Ntrials = 1,
                    num.threads=1)
  mod
}

# Function to fit validation model with INLA
#  species: species name
#  dataname: data name, string
#  Envnames: vector of names of environmental variables
#  mcmc.coda: mcmc output, as coda list
#  maxent.mod: maxEnt model
#  nclust: number of clusters to use
FitValidationModelINLA <- function(species, dataname, Envnames, 
                                   mcmc.coda, maxent.mod, nclust=2) {
  require(disdat)
  require(INLA)
  require(future.apply)
  # Function to fit calibration model to one MCMC draw
  FitCalib <- function(mcmcpars, env, pres) {
    LP <- env%*%mcmcpars
    LP.c <- LP - mean(LP)
    dat <- data.frame(Pres = pres, LP = LP.c)
    dat$Pred <- dat$LP - mean(dat$LP)
    mod <- INLA::inla(Pres ~ Pred + offset(Pred), data=dat, family = "binomial", Ntrials = 1,
                      num.threads=1)
    mod
  }
  
  # Format the data
  if(dataname%in%c("NSW", "AWT")) {
    Potmp <- disPo(dataname)
    grp <- Potmp$group[Potmp$spid==species][1]
    PA <- disPa(dataname, group = grp)[,c("siteid", species)]
    PAEnv <- disEnv(dataname, group = grp)
  } else {
    PA <- disPa(dataname)[,c("siteid", species)]
    PAEnv <- disEnv(dataname)
  }
  
  # scale covariates, using PO data
  EnvScales <- apply(disPo(dataname)[,Envnames], 2, function(x) c(mean=mean(x), sd=sd(x)))
  PAEnv[,Envnames] <- scale(PAEnv[,Envnames], center = EnvScales["mean",], scale=EnvScales["sd",])
  SpData.val <- merge(PA, PAEnv, by="siteid")
  Env <- as.matrix(PAEnv[,Envnames])
  Pres <- SpData.val[,species]
  
  # Validate MaxEnt model
  mxDat <- data.frame(
    Pred =  predict(maxent.mod, SpData.val, type="link"),
    Pres = Pres
  )
  
  # need formula
  f <- formula(paste0("~ ", paste(maxent.mod$beta@Dimnames[[1]], collapse=" + "), "-1"))
  preddat <- model.matrix(f, data.frame(SpData.val))
  
  mxDat <-data.frame(
    Pred = preddat[,names(maxent.mod$betas)]%*%maxent.mod$betas,
    #    Pred.ME = predict(maxent.mod, SpData.val, type="link"), # the same up to constant
    Pres = Pres
  )
  
  mxDat$Pred <- mxDat$Pred - mean(mxDat$Pred) # centre predicitons
  MaxNet.val <- inla(Pres ~ Pred + offset(Pred), data=mxDat, Ntrials = 1, 
                     family="binomial", num.threads=nclust)
  
  mcmc <- as.matrix(mcmc.coda)[,grep("beta", colnames(mcmc.coda[[1]]))]
  colnames(mcmc) <- colnames(preddat)
  
  # Validate Bayesian model
  # parallelise with future
  models <- future_sapply(1:nrow(mcmc), function(X, mcmc, env, pres) {
    FitCalib(mcmc[X,], env=env, pres=pres)
  }, mcmc=mcmc, env=preddat, pres=Pres, simplify = FALSE)
  
  # Merge INLA results
  merged <- inla.merge(models)
  
  Quantiles <- simplify2array(lapply(merged$marginals.fixed, 
                                     function(m) inla.qmarginal(c(0.025, 0.5, 0.975), m)))
  rownames(Quantiles) <- c("0.025quant", "0.5quant", "0.975quant")
  
  GetMarg <- function(marg, fn=function(x) x) {
    unlist(lapply(marg, function(m, ff) inla.emarginal(ff, m),
                  ff=fn))
  }
  Ex <- GetMarg(marg=merged$marginals.fixed, fn=function(x) x)
  Ex2 <- GetMarg(marg=merged$marginals.fixed, fn=function(x) x^2)
  
  Summ <- t(rbind(mean = Ex, sd = sqrt(Ex2 - Ex^2), 
                  mode=lapply(merged$marginals.fixed, inla.mmarginal),
                  Quantiles,
                  kld=rep(0, length(merged$marginals.fixed))))
  
  res <- list(MaxNet.summ = summary(MaxNet.val)$fixed,
              BayesMaxEnt = Summ)
  res
}


# Validate model
# Arguments:
#  sp: species name
#  DataName: data name, string
#  Env: vector of names of environmental variables
#  nclust: number of clusters to use
#  small: If TRUE will run 2 short chains, for testing. Defaults to FALSE
#  classes: Classes of MaxEnt features to use. Defaluts to "l", i.e. just linear
#  verbose: Should the function tell us what it has done? Defaults to FALSE
#  fileprefix: Prefix (incl. folders) for in & output. Defaults to "Results/MaxEntRes" for compatability
JustValidate <- function(sp, DataName, Env, nclust=1, small=FALSE, classes="l", 
                         verbose=FALSE, overwrite=FALSE, fileprefix="Results/MaxEntRes") {
  InName <- paste0(fileprefix, "MCMC", DataName, "_", classes, "_", sp, ".RData")
  OutName <- paste0(fileprefix, "INLA", DataName, "_", classes, "_", sp, ".RData")
  if(small) {
    InName <- gsub("MCMC", "MCMCsmall", InName)
    OutName <- gsub("INLA", "INLAsmall", InName)
  }
  if(!file.exists(OutName)|overwrite) {
    if(verbose) message("Started ", sp)
    require(disdat)
    
    load(InName)
    Valid <- FitValidationModelINLA(species=sp, dataname=DataName, 
                                    Envnames=Env,
                                    mcmc.coda=Fit$mcmc, maxent.mod=Fit$maxnet, nclust=nclust)
    save(Valid, file=OutName)
    if(verbose) message("Done ", sp)
  } else {
    if(verbose) message("Already Done ", sp)
  }
}


# Function to write bash script that will run separate fits for each species
WriteBashScript <- function(dataname, filename, classes="l") {
  cat("#!/bin/bash\n", file=filename)
  Species <- unique(disPo(dataname)$spid)
  
  sapply(Species, function(sp, datname) {
    Command <- paste("Rscript MaxEntBayesScriptSpecies.R --args ", 
                     datname, sp,  classes, "&\n", sep=" ")
    cat(Command, file=filename, append = TRUE)
  }, datname=dataname)
}




















# Fit and validate model
# FitAndValidate <- function(sp, DataName, Env, nclust=1, small=FALSE, classes="l", verbose=FALSE) {
#   FileName <- paste0("Results/MaxEntResINLA", DataName, "_", classes, "_", sp, ".RData")
#   if(small) {
#     NBurn <- 10
#     NIter <- 10
#     NChain <- 2
#     Thin=1
#     FileName <- gsub("INLA", "INLAsmall", FileName)
#   } else {
#     NBurn <- 1.5e4
#     NIter <- 1e3
#     NChain <- 5
#     Thin <- 5
#   }
#   if(!file.exists(FileName)) {
#     if(verbose) message("Started ", sp)
#     require(disdat)
#     Data <- disData(DataName)
#     Fit <- FitModelToDisDat(species=sp, data=Data, Envnames=Env, classes=classes,
#                             nburn =NBurn, niter = NIter, nchain=NChain, thin=Thin)
#     if(verbose) cat("Fitted ", sp, "\n")
#     
#     Valid <- FitValidationModelINLA(species=sp, dataname=DataName, 
#                                     Envnames=Env,
#                                     mcmc.coda=Fit$mcmc, maxent.mod=Fit$maxnet, nclust=nclust)
#     save(Fit, Valid, file=FileName)
#     if(verbose) message("Done ", sp)
#   } else {
#     if(verbose) message("Already Done ", sp)
#   }
# }


# Extract posteriors for intercept and slope
GetResults <- function(sp, dataname, classes="l", fileprefix="Results/MaxEntRes") {
  Filemane <- paste0(fileprefix, "INLA", dataname, "_", classes, "_", sp, ".RData")
  if(file.exists(Filemane)) {
    load(Filemane)
    GetStats <- function(summ, prefix="") {
      summ2 <- summ[,c("mean", "sd")]
      stats.c <- c(summ2)
      names(stats.c) <- paste(prefix, 
                              rep(rownames(summ2), times=2),
                              rep(colnames(summ2), each=2), sep=".")
      
      stats.c
    }
    res <- unlist(c(GetStats(Valid$MaxNet.summ, prefix="MaxEnt"), 
                    GetStats(Valid$BayesMaxEnt, prefix="BayesMaxEnt")))
    if(is.list(res)) res <- unlist(res)
    res
  }
}

# Extract/calculate limits for plots
CalcLims <- function(df, modname) {
  mn <- grep(paste0("^", modname, ".mean"), names(df))
  sd <- grep(paste0("^", modname, ".sd"), names(df))
  
  res <- cbind(
    LLim = unlist(df[,mn]) - df[,sd],
    ULim = df[,mn] + df[,sd]
  )
  colnames(res) <- paste(modname, colnames(res), sep=".")
  res
}  


# Function to run the fitting and or validation
# dataname: name of data in disdat
# remove: variables to remove from environmental dataframe
# classes: which classes to use in MaxEnt. Defaults to "l",
# fit: Should the fitting be done?
#  fileprefix: Prefix (incl. folders) for in & output. Defaults to "Results/MaxEntRes" for compatability
# validate: Should the validation be done?
FitAll <- function(dataname, remove, classes="l", fit=TRUE, validate=TRUE, lambda=NULL, 
                   overwrite=FALSE, fileprefix="Results/MaxEntRes") {
  if(!fit & !validate) stop("fit & validate both FALSE, so stopping here")
  Env <- names(disBg(dataname))
  Env <- Env[!(Env%in%remove)]
  Species <- unique(disPo(dataname)$spid)
  
  if(fit) Thing1 <- sapply(Species, JustFit, DataName=dataname, Env=Env, lambda=lambda, 
                           classes=classes, overwrite=overwrite, fileprefix=fileprefix)
  if(validate) Thing2 <- sapply(Species, JustValidate, DataName=dataname, 
                                Env=Env, small=FALSE, nclust = 6, classes="l", 
                                verbose=TRUE, overwrite=overwrite, fileprefix=fileprefix)
}


# Calcualte goodness of fit statistics
CalcFitStats <- function(pres, pred, thresh=NULL) {
  require(Metrics)
  if(is.null(thresh)) thresh <- mean(pred)
  PredOne <- pred > thresh
  
  #  ConfMat <- table(pres, PredOne)
  Sens <- sum(PredOne*pres)/sum(PredOne)
  Spec <- sum((1-PredOne)*(1-pres))/sum(1-PredOne)
  TSS <- Sens + Spec - 1
  AUC <- Metrics::auc(pres, pred)
  c(Sensitivity = Sens, Specificity = Spec, TSS=TSS, AUC=AUC)  
}


# Fit MaxEnt to one species
# sp: code for species
# dataname: name of data in disdat
#  EnvNames: vector of names of environmental variables
# classes: which classes to use in MaxEnt. Defaults to "l",
# verbose: Should the function update how far it has got?
# valid: Should validation statistics be calculated? Default: FALSE, 
# pred: Should predictions be returned? Default: FALSE
# prob: Should the prediction be on the probability scale? Default: FALSE
# link: link function. Default: "logit", 
# Should records of other other species be used for background points? Default:  FALSE

FitMaxEntToSp <- function(sp, dataname, EnvNames, classes, verbose=FALSE, 
                          valid=FALSE, pred=FALSE, link="logit", otherSpBG = FALSE, prob = FALSE) {
  require(glmnet)
  if(verbose) message("Starting ", sp)
  Pres <- disPo(dataname)[disPo(dataname)$spid==sp,c("occ", EnvNames)]
  if(otherSpBG) {
    bgEnv <- disPo(dataname)[disPo(dataname)$spid!=sp,]
    bgEnv$SiteLoc <- paste0(bgEnv$x, bgEnv$y)
    Use <- sapply(unique(bgEnv$SiteLoc), 
                  function(site, dat) which(dat$SiteLoc==site)[1], dat=bgEnv)
    bgEnv <- bgEnv[Use,]
    bgEnv$occ <- 0 # otherwise everything is a presence...
  } else {
    bgEnv <- disBg(dataname)
  }
  
  dat <- rbind(Pres, bgEnv[,c("occ", EnvNames)])
  if(dataname=="CAN") dat$ontveg <- factor(dat$ontveg)
  if(dataname=="NSW") dat$vegsys <- factor(dat$vegsys)
  if(dataname=="NZ") {
    dat$age <- factor(dat$age)
    dat$toxicats <- factor(dat$toxicats)
  }
  if(dataname=="SWI") dat$calc <- factor(dat$calc)
  
  EnvNamesNoFactor <- EnvNames[!sapply(dat[,EnvNames], is.factor)]
  EnvScales <- apply(dat[,EnvNamesNoFactor], 2, function(x) c(mean=mean(x), sd=sd(x)))
  dat[,EnvNamesNoFactor] <- scale(dat[,EnvNamesNoFactor], 
                                  center = EnvScales["mean",], 
                                  scale=EnvScales["sd",])
  
  #  Fit MaxEnt model
  MaxNet.mod <- try(maxnet(p=dat$occ, data=dat[,EnvNames],
                           f=maxnet.formula(p=dat$occ, data=dat[,EnvNames],
                                            classes=classes))) 
  if(verbose) message("Maxent for ", sp, " done")
  
  # Validate on PA data
  
  if(dataname%in%c("NSW", "AWT")) {
    Potmp <- disPo(dataname)
    grp <- Potmp$group[Potmp$spid==sp][1]
    PA <- disPa(dataname, group = grp)[,c("siteid", sp)]
    PAEnv <- disEnv(dataname, group = grp)
  } else {
    PA <- disPa(dataname)[,c("siteid", sp)]
    PAEnv <- disEnv(dataname)
  }
  PA[,sp] <- factor(PA[,sp], levels=c(0,1))
  
  if(dataname=="CAN") PAEnv$ontveg <- factor(PAEnv$ontveg)
  if(dataname=="NSW") PAEnv$vegsys <- factor(PAEnv$vegsys)
  if(dataname=="NZ") {
    PAEnv$age <- factor(PAEnv$age)
    PAEnv$toxicats <- factor(PAEnv$toxicats)
  }
  if(dataname=="SWI") PAEnv$calc <- factor(PAEnv$calc)
  
  PAEnv[,EnvNamesNoFactor] <- scale(PAEnv[,EnvNamesNoFactor], 
                                    #                                    center = EnvScales["mean",], 
                                    center = TRUE, 
                                    scale=EnvScales["sd",])
  Validdata <- merge(PA, PAEnv, by="siteid")
  
  # Validate MaxEnt model
  Validdata$Pred <-  predict(MaxNet.mod, Validdata, type="link")
  Validdata$Pred <-  Validdata$Pred - mean(Validdata$Pred)
  f <- formula(paste0(sp, " ~ Pred"))
  
  #  PAmod <- glm(f.glm, data=Validdata, family=binomial(link))
  validmod <- glm(f, data=Validdata, family=binomial(link))
  p.p <- as.numeric(as.character(Validdata[,sp]))
  
  PAmod <- maxnet(p=p.p, data=Validdata[,EnvNames], wt=1, 
                  f=maxnet.formula(p=p.p, 
                                   data=Validdata[,EnvNames],
                                   classes=classes))
  
  if(verbose) message("Validation for ", sp, " done")
  # Combine predictions  
  if(pred | valid) {
    Predicted = data.frame(maxnet = Validdata$Pred,
                           valid = predict(validmod, newdata=Validdata, type = "link"),
                           PA = predict(PAmod, newdata=Validdata[,EnvNames], type = "link")
    )
    if(prob) {
      Predicted.prob = data.frame(maxnet = predict(MaxNet.mod, Validdata, type="logistic"),
                                  valid = predict(validmod, newdata=Validdata, type = "response"),
                                  PA = predict(PAmod, newdata=Validdata[,EnvNames], type = "logistic")
      )
    }
    Predicted <- apply(Predicted, 2, scale, scale=FALSE)
  }
  # Calculate validation Statistics
  if(valid) {
    Valid <- apply(Predicted, 2, function(pred, pres)  {
      CalcFitStats(pres=pres, pred=pred, thresh=NULL)
    }, pres=as.numeric(as.character(Validdata[,sp])))
  }
  res <- list(coefficients = coef(validmod))
  if(pred) {
    res$pred <- Predicted
    if(prob) res$pred.prob <- Predicted.prob
  }
  if(valid) res$valid <- Valid
  res
}



# Just fit & validate maxEnt
JustMaxEnt <- function(dataname, remove, classes="l", verbose=FALSE, 
                       link="logit", ...) {
  require(disdat)
  require(maxnet)
  require(future.apply)
  bgEnv <- disBg(dataname)
  EnvNames <- names(bgEnv)[!(names(bgEnv)%in%remove)]
  SpNames <- unique(disPo(dataname)$spid)
  
  #  sp <- SpNames[38]
  
  Coefs <- future_sapply(SpNames, FitMaxEntToSp, dataname=dataname, 
                         EnvNames=EnvNames, classes=classes, verbose=verbose, 
                         future.seed=TRUE, link=link, simplify=FALSE, ...)
  Coefs
}

