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
                          regmult = 1, lambda=1, classes=classes)
  
  # Fit with MaxNet. This will give us the optimum lambda
  MaxNet.mod <- maxnet(p=ToNimble$Data$y, data=as.data.frame(ToNimble$Data$X), 
                       f=maxnet.formula(p=ToNimble$Data$y, data=as.data.frame(ToNimble$Data$X), 
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

  mcmc <- as.matrix(mcmc.coda)[,grep("beta", colnames(mcmc.coda[[1]]))]
  colnames(mcmc) <- Envnames

  # Validate MaxEnt model
  mxDat <- data.frame(
    Pred =  predict(maxent.mod, SpData.val, type="link"),
    Pres = Pres
  )
  mxDat$Pred <- mxDat$Pred - mean(mxDat$Pred)
  MaxNet.val <- inla(Pres ~ Pred + offset(Pred), data=mxDat, Ntrials = 1, 
                     family="binomial", num.threads=nclust)
  
  # Validate Bayesian model
  this_cluster <- parallel::makeCluster(nclust)

  parallel::clusterExport(cl = this_cluster,
                          varlist = c("FitCalib", "mcmc")) # "Env", "Pres"))

  models <- parSapplyLB(this_cluster, 1:nrow(mcmc), function(X, mcmc, env, pres) {
    FitCalib(mcmc[X,], env=env, pres=pres)
  }, mcmc=mcmc, env=Env, pres=Pres, simplify = FALSE)
  
  parallel::stopCluster(this_cluster)
  
# Merge INLA results
  Bayes.val <- do.call(merge, lapply(models, function(l) l))
  
  res <- list(MaxNet.summ = summary(MaxNet.val)$fixed,
              BayesMaxEnt = summary(Bayes.val)$fixed)
  res
}


# Fit and validate model
FitAndValidate <- function(sp, DataName, Env, nclust=1, small=FALSE, classes="l", verbose=FALSE) {
  FileName <- paste0("Results/MaxEntResINLA", DataName, sp, ".RData")
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
    if(verbose) cat("Started ", sp, "\n")
    require(disdat)
    Data <- disData(DataName)
    Fit <- FitModelToDisDat(species=sp, data=Data, Envnames=Env, classes=classes,
                            nburn =NBurn, niter = NIter, nchain=NChain, thin=Thin)
    if(verbose) cat("Fitted ", sp, "\n")
    
    Valid <- FitValidationModelINLA(species=sp, dataname=DataName, 
                                    Envnames=Env,
                                    mcmc.coda=Fit$mcmc, maxent.mod=Fit$maxnet, nclust=nclust)
    save(Fit, Valid, file=FileName)
    if(verbose) cat("Done ", sp, "\n")
  } else {
    if(verbose) cat("Already Done ", sp, "\n")
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
