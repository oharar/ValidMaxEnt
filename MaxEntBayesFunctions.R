# Function to fit Bayesian MaxEnt model to a sopecies in a DisDat dataset
# species: species name/code
# data: DisDat data
# Envnames: vector of names of environmental variables
# nburn: Number of burn-in iterations, defaults to 10
# niter: Number of iterations after burn-in, defaults to 10
# nchain: Number of chains, defaults to 1
FitModelToDisDat <- function(species, data, Envnames, 
                             nburn =10, niter = 10, nchain=1, thin=1, classes="l") {
  SpData <- rbind(data$po[data$po$spid==species,c("occ", Envnames)], data$bg[,c("occ", Envnames)])
  EnvScales <- apply(SpData[,Envnames], 2, function(x) c(mean=mean(x), sd=sd(x))) # scale covariates
  SpData[,Envnames] <- apply(SpData[,Envnames], 2, scale) # scale covariates
  
  ToNimble <- SetUpMaxEnt(p=SpData$occ, data=SpData[,Envnames], 
                          regmult = 1, lambda=1, classes=classes, 
                          addsamplestobackground=FALSE)
  
  # Fit with MaxNet. This will give us the optimum lambda
  MaxNet.mod <- maxnet(p=ToNimble$Data$y, data=SpData[,Envnames], 
                       f=maxnet.formula(p=ToNimble$Data$y, data=SpData[,Envnames], 
                                      classes=classes))
  
  
  ToNimble$Const$MeanLambda <- MaxNet.mod$lambda[200]
  
  output <- FitMaxEnt(maxdat=ToNimble, parallel = TRUE, adaptInterval=nburn,
                      nchains=nchain, nburnin = nburn, niter=nburn+niter, thin=thin)
  list(maxnet = MaxNet.mod, mcmc = output)
}



FitCalib <- function(mcmcpars, env, pres) {
  LP <- env%*%mcmcpars
  LP.c <- LP - mean(LP)
  dat <- data.frame(Pres = pres, LP = LP.c)
  mod <- INLA::inla(Pres ~ LP, data=dat, family = "binomial", Ntrials = 1)
  mod
}


# Function to fit validation model with INLA
FitValidationModelINLA <- function(species, dataname, Envnames, 
                                   mcmc.coda, maxent.mod, nclust=2) {
  require(disdat)
  require(INLA)
  require(parallel)
# Function to fit calibration model to one MCMC draw
  FitCalib <- function(mcmcpars, env, pres) {
    LP <- env%*%mcmcpars
    LP.c <- LP - mean(LP)
    dat <- data.frame(Pres = pres, LP = LP.c)
    dat$LP <- dat$LP - mean(dat$LP)
    mod <- INLA::inla(Pres ~ LP + offset(LP), data=dat, family = "binomial", Ntrials = 1,
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

  mxDat$Pred <- mxDat$Pred - mean(mxDat$Pred)
  MaxNet.val <- inla(Pres ~ Pred + offset(Pred), data=mxDat, Ntrials = 1, 
                     family="binomial", num.threads=nclust)
  
  mcmc <- as.matrix(mcmc.coda)[,grep("beta", colnames(mcmc.coda[[1]]))]
  colnames(mcmc) <- colnames(preddat)
  
  # Validate Bayesian model
  this_cluster <- parallel::makeCluster(nclust)

  parallel::clusterExport(cl = this_cluster,
                          varlist = c("FitCalib", "mcmc", "preddat")) # "Env", "Pres"))

  models <- parSapplyLB(this_cluster, 1:nrow(mcmc), function(X, mcmc, env, pres) {
    FitCalib(mcmc[X,], env=preddat, pres=pres)
  }, mcmc=mcmc, env=Env, pres=Pres, simplify = FALSE)
  
  parallel::stopCluster(this_cluster)
  
# Merge INLA results
#  Bayes.val <- do.call(merge, lapply(models, function(l) l))
  merged <- inla.merge(models)
  
  # sd 0.025quant 0.5quant 0.975quant  
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


# Fit model
JustFit <- function(sp, DataName, Env, small=FALSE, classes="l", verbose=FALSE) {
  FileName <- paste0("Results/MaxEntResMCMC", DataName, "_", classes, "_", sp, ".RData")
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
  if(!file.exists(FileName)) {
    if(verbose) message("Started ", sp)
    require(disdat)
    Data <- disData(DataName)
    Fit <- FitModelToDisDat(species=sp, data=Data, Envnames=Env, classes=classes,
                            nburn =NBurn, niter = NIter, nchain=NChain, thin=Thin)
    if(verbose) message("Fitted ", sp)
    
    save(Fit, file=FileName)
    if(verbose) message("Done ", sp)
  } else {
    if(verbose) message("Already Done ", sp)
  }
}

# Validate model
JustValidate <- function(sp, DataName, Env, nclust=1, small=FALSE, classes="l", verbose=FALSE) {
  InName <- paste0("Results/MaxEntResMCMC", DataName, "_", classes, "_", sp, ".RData")
  OutName <- paste0("Results/MaxEntResINLA", DataName, "_", classes, "_", sp, ".RData")
  if(small) {
    InName <- gsub("MCMC", "MCMCsmall", InName)
    OutName <- gsub("INLA", "INLAsmall", InName)
  }
  if(!file.exists(OutName)) {
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


# Fit and validate model
FitAndValidate <- function(sp, DataName, Env, nclust=1, small=FALSE, classes="l", verbose=FALSE) {
  FileName <- paste0("Results/MaxEntResINLA", DataName, "_", classes, "_", sp, ".RData")
  if(small) {
    NBurn <- 10
    NIter <- 10
    NChain <- 2
    Thin=1
    FileName <- gsub("INLA", "INLAsmall", FileName)
  } else {
    NBurn <- 1.5e4
    NIter <- 1e3
    NChain <- 5
    Thin <- 5
  }
  if(!file.exists(FileName)) {
    if(verbose) message("Started ", sp)
    require(disdat)
    Data <- disData(DataName)
    Fit <- FitModelToDisDat(species=sp, data=Data, Envnames=Env, classes=classes,
                            nburn =NBurn, niter = NIter, nchain=NChain, thin=Thin)
    if(verbose) cat("Fitted ", sp, "\n")
    
    Valid <- FitValidationModelINLA(species=sp, dataname=DataName, 
                                    Envnames=Env,
                                    mcmc.coda=Fit$mcmc, maxent.mod=Fit$maxnet, nclust=nclust)
    save(Fit, Valid, file=FileName)
    if(verbose) message("Done ", sp)
  } else {
    if(verbose) message("Already Done ", sp)
  }
}


# Extract posteriors for intercept and slope
GetResults <- function(sp, dataname) {
  Filemane <- paste0("MaxEntResults/MaxEntResINLA", dataname, sp, ".RData")
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
    res <- c(GetStats(Valid$MaxNet.summ, prefix="MaxEnt"), 
             GetStats(Valid$BayesMaxEnt, prefix="BayesMaxEnt"))
    res
  }
}

# Extract/calculate limits for plots
CalcLims <- function(df, modname) {
  mn <- grep(paste0("^", modname, ".mean"), names(df))
  sd <- grep(paste0("^", modname, ".sd"), names(df))
  
  res <- cbind(
    LLim = df[,mn] - df[,sd],
    ULim = df[,mn] + df[,sd]
  )
  colnames(res) <- paste(modname, colnames(res), sep=".")
  res
}  
