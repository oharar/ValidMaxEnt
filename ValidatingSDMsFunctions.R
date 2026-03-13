
# Plot a species' data
PlotMap <- function(sp) {
  require(sf)
  require(disdat)
  region <- toupper(gsub("[0-9]*", "", sp))
  bg <- disBg(region)
  po <- disPo(region)
  
  if (region %in% c("AWT", "NSW")) {
    grp <- po[po$spid == sp, "group"][1]
    pa <- disPa(region, grp)
  } else {
    pa <- disPa(region)
  }
  
  plot(sf::st_geometry(disBorder(region)), main=sp)
  points(bg$x, bg$y, col="grey80", cex=0.2, pch=16)
  points(pa$x[pa[,sp]==0], pa$y[pa[,sp]==0], col="pink", cex=0.3, pch=16)
  points(pa$x[pa[,sp]==1], pa$y[pa[,sp]==1], col="red", cex=0.5, pch=16)
  points(po$x[po$spid==sp], po$y[po$spid==sp], col="blue", cex=0.5, pch=16)
}







# Calcualte goodness of fit statistics
CalcFitStats <- function(pres, pred, thresh=NULL) {
  require(Metrics)
  if(is.null(thresh)) thresh <- mean(pred)
  PredOne <- pred > thresh
  
  #  ConfMat <- table(pres, PredOne)
  Sens <- sum(PredOne*pres)/sum(PredOne + 1e-6)
  Spec <- sum((1-PredOne)*(1-pres))/sum(1-PredOne + 1e-6)
  TSS <- Sens + Spec - 1
  AUC <- Metrics::auc(pres, pred)
  c(Sensitivity = Sens, Specificity = Spec, TSS=TSS, AUC=AUC)  
}


# Fit MaxEnt to one species
# sp: code for species
# region: name of data in disdat
#  EnvNames: vector of names of environmental variables
# classes: which classes to use in MaxEnt. Defaults to "l",
# verbose: Should the function update how far it has got?
# valid: Should validation statistics be calculated? Default: FALSE, 
# pred: Should predictions be returned? Default: FALSE
# prob: Should the prediction be on the probability scale? Default: FALSE
# link: link function. Default: "logit", 
# Should records of other other species be used for background points? Default:  FALSE

FitMaxEntToSp <- function(sp, classes, verbose=FALSE, link="logit", 
                          removenames=c("siteid", "spid", "x", "y", "occ", "group"), 
                          valid=FALSE, pred=FALSE, otherSpBG = FALSE, prob = FALSE, 
                          savemodels = FALSE) {
  require(glmnet)
  if(verbose) message("Starting ", sp)
  region <- toupper(gsub("[0-9]*", "", sp))
  
  # Get presence-only data
  POdata <- GetPOdata(sp=sp, otherSpBG=otherSpBG, scale=TRUE)
  EnvNames <- names(POdata)[!names(POdata)%in%removenames]
  #  Fit MaxEnt model
  MaxNet.mod <- FitMaxEnt(dat=POdata, EnvNames=EnvNames, RespName="occ", 
                          PA=FALSE, classes=classes)
  
  if(verbose) message("Maxent for ", sp, " done")
  
  # Get PA data and centre
  PAdata <- GetPAdata(sp=sp, scale=FALSE)
  PAdata$p.p <- as.numeric(as.character(PAdata$PresAbs)) # convert PA to numeric from logical
  
  # better to pass these into GetPAdata()?
  if(!is.null(attr(POdata, "envmeans"))) {
    PAdata[,names(attr(POdata, "envmeans"))] <- sweep(PAdata[,names(attr(POdata, "envmeans"))], 
                                                      2, attr(POdata, "envmeans"), "-")
    PAdata[,names(attr(POdata, "envsds"))] <- sweep(PAdata[,names(attr(POdata, "envsds"))], 
                                                    2, attr(POdata, "envsds"), "/")
  }
  
  # Validate MaxEnt model, and calculate predictions. 
  # If the model dodn't work, create lots of NAs
  if(!is.null(MaxNet.mod)) {
    PAdata$MxPred <-  predict(MaxNet.mod, PAdata, type="link")
    PAdata$MxPred <-  PAdata$MxPred - mean(PAdata$MxPred)
    validmod <- glm(PresAbs~MxPred, data=PAdata, family=binomial(link))
    PAdata$valid <- predict(validmod, newdata=PAdata, type = "link")
    
    if(prob) {
      PAdata$MxProb <- predict(MaxNet.mod, PAdata, type="logistic")
      PAdata$validProb <- predict(validmod, newdata=PAdata, type = "response")
    }
    # Create results: more terms get added later
    res <- list(coefficients = coef(validmod),
                alpha = MaxNet.mod$alpha, 
                confint = confint(validmod))
    if(savemodels) {
      res$MaxEnt <- MaxNet.mod
      res$valid <- validmod
    }
    
  } else {
    PAdata$MxPred <- PAdata$valid <- PAdata$validProb <- PAdata$MxProb <- rep(NA,nrow(PAdata))
    res <- list(coefficients = NA, alpha = NA, confint = NA)
  }
  
  # Fit MaxEnt model to PA data as PA
  
  PAmod <- FitMaxEnt(dat=PAdata, EnvNames=EnvNames, RespName="p.p", 
                     PA=TRUE, classes=classes)
  if(!is.null(PAmod)) {
    PAdata$PApred <- predict(PAmod, newdata=PAdata[,EnvNames], type = "link")
    if(prob) {
      PAdata$PAprob <- predict(PAmod, newdata=PAdata[,EnvNames], type = "logistic")
    }
    if(savemodels) res$PAMaxEnt <- PAmod
    
  } else {
    PAdata$PApred <- rep(NA,nrow(PAdata))
  }
  
  if(verbose) message("Validation for ", sp, " done")
  # Combine predictions  
  if(pred | valid) {
    if(!is.null(PAmod)) {
      Predicted <- data.frame(maxnet = PAdata$MxPred,
                              valid = PAdata$valid,
                              PA = PAdata$PApred
      )  
      Predicted <- apply(Predicted, 2, scale, scale=FALSE) # Mean centre
      res$pred <- Predicted
      
    } else {
      Predicted <- NULL
    }
    
    # Return probabilities
    if(prob) {
      Predicted.prob <- data.frame(maxnet = PAdata$MxProb,
                                   valid = PAdata$validProb,
                                   PA = c(PAdata$PAprob))
      res$pred.prob <- Predicted.prob
    }
    
  }
  
  # Calculate validation Statistics
  if(valid & !is.null(Predicted)) {
    Valid <- apply(Predicted, 2, function(pred, pres)  {
      res$valid <- CalcFitStats(pres=pres, pred=pred, thresh=NULL)
    }, pres=PAdata$p.p)
    res$Valid <- Valid
  }
  return(res)
}

# Get data and fit models. 
# sp: species
# classes: classes of MaxEnt model. Defaults to "l"
# verbose: Should progress be reported? Defaults to FALSE
# sim: Should simulated data be used? Defaults to FALSE,
# future.seed: Should this be future-proved for aprallelisation?
# link: link function. Default: logit=link, 
#  ...: other arguments thrown around
GetRes <- function(sp, classes="l", verbose=verbose, sim=FALSE, 
                   future.seed=TRUE, link="logit", onlySumm = FALSE, ...) {
  if(sim) {
    Data <- SimStatsFromData(sp) #, ...)
  } else {
    Data <- GetData(sp, ...)
  }
  FM <- FitModels(Data, classes=classes, verbose=verbose) #, ...)  

  if(onlySumm) {
    res <- c(corModels = ifelse(is.null(FM$pred), NA, cor(FM$pred[,"maxnet"]), c(FM$pred[,"PA"])),
             AUC = ifelse(is.null(FM$valid), NA, FM$valid["AUC", "maxnet"]),
             TSS = ifelse(is.null(FM$pred), NA, FM$valid["TSS", "maxnet"]),
             beta = ifelse(is.null(FM$coefficients), NA, FM$coefficients["MxPred"]))
    return(res)
  } else {
    return(FM)
  }
}

# Just fit & validate maxEnt
JustMaxEnt <- function(region, remove=c("siteid", "spid", "x", "y", "occ", "group"), 
                       classes="l", verbose=FALSE, POthreshold = 20, 
                       link="logit", ...) {
  require(disdat)
  require(maxnet)
  require(future.apply)
  bgEnv <- disBg(region)
  EnvNames <- names(bgEnv)[!(names(bgEnv)%in%remove)]
  Spp <- table(disPo(region)$spid)
  SpNames <- names(Spp)[Spp>=POthreshold]
  
  #  sp <- SpNames[38]
  #  Coefs <- future_sapply(SpNames, FitMaxEntToSp, 
  Coefs <- future_sapply(SpNames, GetRes, 
                         classes=classes, verbose=verbose, 
                         future.seed=TRUE, link=link, simplify=FALSE, ...)
  Coefs
}

# Just fit & validate maxEnt
JustSimMaxEnt <- function(region, remove=c("siteid", "spid", "x", "y", "occ", "group"), 
                          classes="l", verbose=FALSE, link="logit", POthreshold = 20, ...) {
  require(disdat)
  require(maxnet)
  require(future.apply)
  bgEnv <- disBg(region)
  EnvNames <- names(bgEnv)[!(names(bgEnv)%in%remove)]
  Spp <- table(disPo(region)$spid)
  SpNames <- names(Spp)[Spp>=POthreshold]
  
  #  sp <- SpNames[38]
  #  Coefs <- future_sapply(SpNames, FitMaxEntToSp, 
  Coefs <- future_sapply(SpNames, GetRes, sim=TRUE, 
                         classes=classes, verbose=verbose, 
                         future.seed=TRUE, link=link, simplify=FALSE, ...)
  Coefs
}


# Plot regression coefficients
PlotCoefs <- function(nm, lst, AddPoints=FALSE, AddCIs=TRUE) {
  Coefs.l <- lapply(lst[[nm]], function(l) l$coefficients)
  Coefs <- t(list2DF(Coefs.l)); 
  colnames(Coefs) <- names(lst[[nm]][[1]]$coefficients)
  
  if(AddCIs) {
    CIs.l <- lapply(lst[[nm]], function(l) c(l$confint))
    CIs <- t(list2DF(CIs.l))
    colnames(CIs) <- paste(rep(rownames(lst[[nm]][[1]]$confint), 2), 
                           rep(colnames(lst[[nm]][[1]]$confint), each=2), sep=":")
    PredRange <- range(c(0,1, range(CIs[,grep("Pred", colnames(CIs))])))
    IntRange <- range(CIs[,grep("(Intercept)", colnames(CIs))])
  } else {
    PredRange <- range(c(0,1, range(Coefs[,grep("Pred", colnames(Coefs))])))
    IntRange <- range(Coefs[,"(Intercept)"])
  }
  
  plot(Coefs, xlim=IntRange, ylim=PredRange, type="n", xlab="", ylab="", main=nm) 
  rect(-100, 0, 100, 1, col="mistyrose", border=NA)
  rect(-100, 0, 100, 0.7, col="pink", border=NA)
  if(AddPoints) points(Coefs)
  if(AddCIs) {
    segments(CIs[,"(Intercept):2.5 %"], Coefs[, grep("Pred", colnames(Coefs))], 
             CIs[,"(Intercept):97.5 %"], Coefs[, grep("Pred", colnames(Coefs))])
    segments(Coefs[,"(Intercept)"], CIs[, grep("Pred:2.5 %", colnames(CIs))], 
             Coefs[,"(Intercept)"], CIs[, grep("Pred:97.5 %", colnames(CIs))])
  }
  box()
}




# Utility function to convert covariates to factors
ConvertFactors <- function(dat, region) {
  if(region=="CAN") dat$ontveg <- factor(dat$ontveg)
  if(region=="NSW") {
    # convert pine forest to dry open forest, to avoid errors
    # Pine forest is "exotic" but also rare. For some species it isn't in the PA data
    # because of course a standardised data set should have these annoying problems
    dat$vegsys[dat$vegsys==8] <- 3
    dat$vegsys <- factor(dat$vegsys, levels=1:9)
  }
  if(region=="NZ") {
    dat$age <- factor(dat$age)
    # Classes 2 & 3 are both rare, and according to the documentation, class 3 does not exist.
    dat$toxicats[dat$toxicats==3] <- 2
    dat$toxicats <- factor(dat$toxicats)
  }
  if(region=="SWI") dat$calc <- factor(dat$calc)
  dat
}


# Function to get presence/absence data for one species
#  returns data frame with: 
#  siteid: site
#  PresAbs: 0: absence, 1: presence
#  group: optional
#  x, y: location(?)
#  then environmental covariates

GetPAdata <- function(sp, removenames = c("siteid", "spid", "x", "y", "occ", "group", "geometry"), 
                      scale=TRUE, sf=FALSE) {
  region <- toupper(gsub("[0-9]*", "", sp))
  
  if(region%in%c("NSW", "AWT")) {
    Potmp <- disPo(region)
    grp <- Potmp$group[Potmp$spid==sp][1]
    PA <- disPa(region, group = grp)[,c("siteid", sp)]
    PAEnv <- disEnv(region, group = grp)
  } else {
    PA <- disPa(region)[,c("siteid", sp)]
    PAEnv <- disEnv(region)
  }
  PAEnv <- ConvertFactors(dat=PAEnv, region=region)
  PA[,sp] <- factor(PA[,sp], levels=c(0,1))
  names(PA)[names(PA)==sp] <- "PresAbs" # change name for consistency
  
  PAdata <- merge(PA, PAEnv, by="siteid")
  
  if(sf) {
    PAdata <- sf::st_as_sf(PAdata, coords = c("x", "y"), crs = disCRS(region)) # subset the species
  }
  
  if(scale) {
    EnvNamesNoFactor <- names(PAdata)[!names(PAdata)%in%removenames & !sapply(PAdata, is.factor)]
    means <- apply(PAdata[,EnvNamesNoFactor, drop=TRUE], 2, mean)
    sds <- apply(PAdata[,EnvNamesNoFactor, drop=TRUE], 2, sd)
    PAdata[,EnvNamesNoFactor] <- scale(PAdata[,EnvNamesNoFactor, drop=TRUE])
    attr(PAdata, "envmeans") <- means
    attr(PAdata, "envsds") <- sds
  }
  
  PAdata
}


# Get PO data
# Returns data frame with "occ, a 0/1 integeter: 0=background, 1=presence, 
#  and columns of the environmental coavriates, converted to factors when appropriate.
GetPOdata <- function(sp, removenames = c("siteid", "spid", "x", "y", "occ", "group", "geometry"), 
                      otherSpBG=FALSE, scale=TRUE, sf=FALSE) {
  region <- toupper(gsub("[0-9]*", "", sp))
  
  if(otherSpBG) {
    bgEnv <- disPo(region)[disPo(region)$spid!=sp,]
    bgEnv$SiteLoc <- paste0(bgEnv$x, bgEnv$y)
    Use <- sapply(unique(bgEnv$SiteLoc), 
                  function(site, dat) which(dat$SiteLoc==site)[1], dat=bgEnv)
    bgEnv <- bgEnv[Use,]
    bgEnv$occ <- 0 # otherwise everything is a presence...
    bgEnv$SiteLoc <- NULL # remove to make consistent
  } else {
    bgEnv <- disBg(region)
  }
  EnvNames <- names(bgEnv)[!(names(bgEnv)%in%removenames)]
  Pres <- disPo(region)[disPo(region)$spid==sp,c("occ", "x", "y", EnvNames)]
  
  dat <- rbind(Pres, bgEnv[,c("occ", "x", "y", EnvNames)])
  dat <- ConvertFactors(dat, region)
  
  if(sf) {
    dat <- sf::st_as_sf(dat, coords = c("x", "y"), crs = disCRS(region)) # subset the species
  }
  if(scale) {
    EnvNamesNoFactor <- names(dat)[!names(dat)%in%removenames & !sapply(dat, is.factor)]
    means <- apply(dat[,EnvNamesNoFactor, drop=TRUE], 2, mean)
    sds <- apply(dat[,EnvNamesNoFactor, drop=TRUE], 2, sd)
    dat[,EnvNamesNoFactor] <- scale(dat[,EnvNamesNoFactor, drop=TRUE])
    attr(dat, "envmeans") <- means
    attr(dat, "envsds") <- sds
  }
  
  dat
}

# Get PO & PA data from disdat, and scale covariates
# sp: specoes
# otherSpBG: Should other species be used as background? Defaults to FALSE
# removenames: variable names that are not environmental. 
#    Default: c("siteid", "spid", "x", "y", "occ", "group")
GetData <- function(sp, otherSpBG = FALSE, scale=TRUE, 
                    removenames=c("siteid", "spid", "x", "y", "occ", "group"), ...) {
  region <- toupper(gsub("[0-9]*", "", sp))
  
  # Get presence-only data
  POdata <- GetPOdata(sp=sp, otherSpBG=otherSpBG, scale=scale)
  EnvNames <- names(POdata)[!names(POdata)%in%removenames]
  
  # Get PA data and centre
  PAdata <- GetPAdata(sp=sp, scale=FALSE)
  PAdata$PresAbsNum <- as.numeric(as.character(PAdata$PresAbs)) # convert PA to numeric from logical
  
  # better to pass these into GetPAdata()?
  if(!is.null(attr(POdata, "envmeans"))) {
    PAdata[,names(attr(POdata, "envmeans"))] <- sweep(PAdata[,names(attr(POdata, "envmeans"))], 
                                                      2, attr(POdata, "envmeans"), "-")
    PAdata[,names(attr(POdata, "envsds"))] <- sweep(PAdata[,names(attr(POdata, "envsds"))], 
                                                    2, attr(POdata, "envsds"), "/")
  }
  res <- list(
    PO = POdata,
    PA = PAdata,
    EnvNames = EnvNames,
    species = sp
  )
  res
}




# Simulate PA or PO data from log weights
# lnWt: log weight (for PO) or probability (for PA)
# PA: should the data be PA or PA? If PA, we assume lnWt is a logit probability
# sigma: standard deviation of overdispersion. Default to NULL, i.e. no overdispersion

SimDataFromWeights <- function(lnWt, PA=FALSE, N=5, sigma=NULL) {
  if(!is.null(sigma)) {
    if(sigma>0) lnWt <- rnorm(length(lnWt), lnWt, sigma)
  }
  if(PA) {
    pr <- 1/(1+exp(-lnWt))
    res <- rbinom(length(lnWt), 1, pr)
  } else {
    Wt <- pmin(1e5,exp(lnWt - mean(lnWt))) # adjust for infinite weights
    SimPO <- sample.int(n=length(Wt), size=N, prob=Wt/sum(Wt), replace=FALSE)
    res  <- as.numeric((1:length(Wt))%in%SimPO)
  }
  res
}



# Function to simulate new data (PA and PO) from model fitted to real data.
# sigma: if not null, add overdispersion (i.e. extra error). Can be length 1 or 2. 
#   If length 2, first element is PO, second is PA
SimFromData <- function(species, sigma=NULL, newdata=NULL, PA=TRUE, 
                        nsim=2, verbose=FALSE, 
                        removenames=c("siteid", "PresAbs", "spid", "x", "y", "occ", "group")) {
  region <- toupper(gsub("[0-9]*", "", species))
  
  require(maxnet)
  if(!is.null(sigma)) {
    if(length(sigma)==1) sigma <- rep(sigma,2)
  } 
  # Get PA data to fit a model to
  Data <- GetData(sp=species, scale=FALSE)
  N <- sum(Data$PO$occ)

  # fit model with linear features to data. If PO, uses (almost) infinite weights.
  YY <- ifelse(PA, "PresAbs ~", "occ ~")
  f <- formula(paste(YY, paste(Data$EnvNames, collapse = " + ")))
  
  if(!PA) {
    Wt <- 1+ (Data$PO$occ * 99)
    truemod <- glm(f, data=Data$PO, family=binomial("logit"), weights=Wt)
  } else {
    Wt <- rep(1, nrow(Data$PA))
    truemod <- glm(f, data=Data$PA, family=binomial("logit"), weights=Wt)
  }
  
  # Simulate data
  # if no new data, simulate on original PA and PO data
  if(is.null(newdata)) {
    Data$PA$eta <- predict(truemod, newdata = Data$PA, type = "link")
    Data$PO$eta <- predict(truemod, newdata = Data$PO, type = "link")
    
    SimPODat <- replicate(nsim, SimDataFromWeights(lnWt=Data$PO$eta, PA=FALSE, N=N, sigma=sigma[1]))
    SimPADat <- replicate(nsim, SimDataFromWeights(lnWt=Data$PA$eta, PA=TRUE, sigma=sigma[2]))
  } else {
    newdata$eta <- predict(truemod, newdata = newdata, type = "link")
    
    SimPODat <- replicate(nsim, SimDataFromWeights(lnWt=newdata$eta, PA=FALSE, N=N, sigma=sigma[1]))
    SimPADat <- replicate(nsim, SimDataFromWeights(lnWt=newdata$eta, PA=TRUE, sigma=sigma[2]))
    
  }
  
  # Simulate data
  
  res <- list(PA=PA, simPA=SimPADat, simPO=SimPODat)
  if(verbose) {
    res$truemod <- truemod
    res$newdata <- newdata
    res$PA <- Data$PA
    res$PO <- Data$PO
    res$EnvNames <- Data$EnvNames
  }
  res
}

# Simulate data from a model
SimDataFromModel <- function(mod, newdata, N, sigma=NULL) {
  if(!is.null(sigma)) {
    if(length(sigma)==1) sigma <- rep(sigma,2)
  }

  if("PA"%in%names(newdata)) {
    eta <- predict(mod, newdata = newdata$PA, type = "link")
    res <- list(PA = newdata$PA)
  } else {
    res <- list(PA = newdata)
    eta <- predict(mod, newdata = newdata, type = "link")
  }
  res$PA$PresAbsNum <- SimDataFromWeights(lnWt=eta, PA=TRUE, sigma=sigma[2])
    
  if("PO"%in%names(newdata)) {
    etaO <- predict(mod, newdata = newdata$PO, type = "link")
    res$PO <- newdata$PO
  } else {
    res$PA <- newdata
    etaO <- predict(mod, newdata = newdata, type = "link")
  }
  res$PO$occ <- SimDataFromWeights(lnWt=etaO, PA=FALSE, N=N, sigma=sigma[1])
  
  res
}


# Function to simulate new data (PA and PO) from model fitted to real data.
# sigma: if not null, add overdispersion (i.e. extra error). Can be length 1 or 2. 
#   If length 2, first element is PO, second is PA
SimStatsFromData <- function(species, sigma=NULL, newdata=NULL, PA=TRUE, ...) {
  region <- toupper(gsub("[0-9]*", "", species))
  
  require(maxnet)
  if(!is.null(sigma)) {
    if(length(sigma)==1) sigma <- rep(sigma,2)
  } 
  # Get PA data to fit a model to
  Data <- GetData(sp=species, scale=FALSE)
  N <- sum(Data$PO$occ)
  if(is.null(newdata)) newdata <- Data
  
  # fit model with linear features to data. If PO, uses (almost) infinite weights.
  YY <- ifelse(PA, "PresAbs ~", "occ ~")
  f <- formula(paste(YY, paste(Data$EnvNames, collapse = " + ")))
  
  if(!PA) {
    Wt <- 1+ (Data$PO$occ * 99)
    truemod <- glm(f, data=Data$PO, family=binomial("logit"), weights=Wt)
  } else {
    Wt <- rep(1, nrow(Data$PA))
    truemod <- glm(f, data=Data$PA, family=binomial("logit"), weights=Wt)
  }
  
  # Simulate data
  # if no new data, simulate on original PA and PO data
  res <- SimDataFromModel(mod=truemod, newdata=newdata, N=N, sigma=sigma)
  res$EnvNames <- Data$EnvNames
  res$species <- species
  
  res
}




# Function to simulate correlations for a species by fitting a GLM to the PA data and 
#  then simulating PA and PO data (PO by simulation on the BG points). 
#  MaxEnt & GLMs are then fitted to the simulated data, and the correlations in 
#  their predictions on the PA data are calculated
#  species <- "swi01"; region <- "SWI"

SimCorrReg <- function(species, nsim=5, stats="corr", 
                       removenames=c("siteid", "spid", "x", "y", "occ", "group"), 
                       ...) {
  require(maxnet)
  region <- toupper(gsub("[0-9]*", "", species))
  
  # Get PA data to fit a model to
  
  Data <- GetData(species, scale = FALSE)
  # PODat <- GetPOdata(species, scale = TRUE)
  # EnvNames <- names(PODat)[!names(PODat)%in%removenames]
  Sim <- SimFromData(species, sigma=sigma, newdata=NULL, nsim=nsim)
  
  Corr <- sapply(1:nsim, GetCorrCoef, PAd=PADat, POd=PODat, sims=Sim, envnames=EnvNames)
  
  Corr
}




# Function to simulate data with calculated weights and probabilities, 
#   and calculate the correlations in the fitted models

# SimCorr <- function(lnWt, N, BG, PA, ENames, stats=c("corr", "reg"), sigma=NULL) {
#   if(!any(stats%in%c("corr", "reg"))) stop("stats must be one or both of corr and reg")
#   if(!is.null(sigma)) lnWt <- rnorm(length(lnWt), lnWt, sigma)
#   Wt <- exp(lnWt - mean(lnWt))
#   SimPO <- sample.int(n=length(Wt), size=N, prob=Wt/sum(Wt))
#   SimOcc  <- as.numeric((1:length(Wt))%in%SimPO)
#   # Fit MaxEnt to PO
#   mod <- maxnet(p=SimOcc, data=BG,
#                        f=maxnet.formula(p=SimOcc, data=BG,
#                                         classes="l"))
#   PA$MxPred <-  predict(mod, PA, type="link")
# 
#   # Fit model to simulated PA data
#   PA$SimPA <- rbinom(nrow(PA), 1, PA$PredProb.true)==1
#   res <- NULL
#   if("corr"%in%stats) {
#     f2 <- formula(paste("SimPA", "~", paste(ENames, collapse = " + ")))
#     GLM.mod <- glm(f2, data=PA, family="binomial")
#     GLMpred <- predict(GLM.mod)
#     res$corr <- cor(GLMpred, MXPred)
#   }
#   if("reg"%in%stats) {
#     GLM.reg <- glm(SimPA ~ MxPred, data=PA, family="binomial")
#     res$reg <- coef(GLM.reg)
#   }
#   res
# }





#  Fit MaxEnt model

#  dat: data
# EnvNames: environmenal variables
# RespName: name of response, defaults tocc, 
# PA: should the data be treated as presence/absence? Defualts to FALSE, so PO
# classes: which classes of model should be sued? Default: "l", i.e. linear

# classes <- "lqpth"
# trydata <- GetPOdata(sp="nz01", region="NZ")
# 
# thing <- FitMxEnt(trydata, RespName = "occ", EnvNames =names(trydata)[!names(trydata)%in% removenames],
#          classes="lqpth")

# dat=PAd; EnvNames=envnames; RespName="simpres"; PA=TRUE; classes="l"
FitMaxEnt <- function(dat, EnvNames, RespName="occ", PA=FALSE, classes="l") {
  
  EnvNamesNotInData <- EnvNames[!sapply(EnvNames, function(EN, dtn) EN%in%dtn, dtn=names(dat))]
  if(length(EnvNamesNotInData)>0) {
    stop(paste0("Variables ", paste(EnvNamesNotInData, collapse=", "), " not in dat"))
  }
  
  wt <- ifelse(PA, 1, 100)
  TryClasses <- sapply(nchar(classes):1, function(wh, str) substr(str, 1, wh), 
                       str=classes)
  
  for(cl in TryClasses) {
    mod <- tryCatch(maxnet(p=dat[,RespName, drop=TRUE], data=dat[,EnvNames, drop=TRUE], wt=wt, 
                           f=maxnet.formula(p=dat[,RespName, drop=TRUE], data=dat[,EnvNames, drop=TRUE],
                                            classes=cl)), 
                    error = function(msg){ return(NULL) })
    if(is.null(mod)) { 
      warning(paste0("Model not converged with ", cl, 
                     " classes. Trying without ", substr(cl, nchar(cl), nchar(cl))))
    } else {
      break
    } 
  }
  mod
}

# Function to get regression coefficient from simulation.
GetRegCoef <- function(ns, PAd, POd, sims, envnames) {
  PAd <- cbind(PAd, simpres=sims$simPA[,ns])
  POd <- cbind(POd, simpres=sims$simPO[,ns])
  
  if(!is.null(attr(POd, "envmeans"))) {
    PAd[,names(attr(POd, "envmeans"))] <- sweep(PAd[,names(attr(POd, "envmeans"))], 
                                                2, attr(POd, "envmeans"), "-")
    PAd[,names(attr(POdata, "envsds"))] <- sweep(PAd[,names(attr(POd, "envsds"))], 
                                                 2, attr(POd, "envsds"), "/")
  }
  
  mod <- FitMaxEnt(dat=POd, EnvNames=envnames, RespName="simpres", PA=FALSE, classes="l")
  PAd$MxPred <-  predict(mod, PAd, type="link")
  PAd$MxPred <-  PAd$MxPred - mean(PAd$MxPred)
  validmod <- glm(simpres~MxPred, data=PAd, family=binomial(logit))
  coef(validmod)[2]
}

# Function to estimate coefficients from simulate PO and PA data based on a species
SimCoefs <- function(species, sigma, nsim, 
                     removenames=c("siteid", "spid", "x", "y", "occ", "group", "simpres")) {
  
  PADat <- GetPAdata(species, scale = FALSE)
  PODat <- GetPOdata(species, scale = TRUE)
  EnvNames <- names(PODat)[!names(PODat)%in%removenames]
  Sim <- SimFromData(species, sigma=sigma, newdata=NULL, nsim=nsim)
  
  Coefs <-  sapply(1:nsim, GetRegCoef, PAd=PADat, POd=PODat, sims=Sim, envnames=EnvNames)
  
  Coefs
  
}

# Get corelation coefficient between simulated PA & PO data
# ns <- 1; PAd=PADat; POd=PODat; sims=Sim; envnames=EnvNames
GetCorrCoef <- function(ns, PAd, POd, sims, envnames) {
  PAd <- cbind(PAd, simpres=sims$simPA[,ns])
  POd <- cbind(POd, simpres=sims$simPO[,ns])
  
  if(!is.null(attr(POd, "envmeans"))) {
    PAd[,names(attr(POd, "envmeans"))] <- sweep(PAd[,names(attr(POd, "envmeans"))], 
                                                2, attr(POd, "envmeans"), "-")
    PAd[,names(attr(POdata, "envsds"))] <- sweep(PAd[,names(attr(POd, "envsds"))], 
                                                 2, attr(POd, "envsds"), "/")
  }
  
  POmod <- FitMaxEnt(dat=POd, EnvNames=envnames, RespName="simpres", PA=FALSE, classes="l")
  PAd$MxPred <-  predict(POmod, PAd, type="link")
  PAmod <- FitMaxEnt(dat=PAd, EnvNames=envnames, RespName="simpres", PA=TRUE, classes="l")
  PAd$Pred <-  predict(PAmod, PAd, type="link")
  cor(PAd$MxPred, PAd$Pred)
}



# Simulate data from model of PO data, and calculate either correlation between 
#  predictions or regression slope for model predicting PA data from PO model
# species: speciees to use
# nsim: number of simulations. Default: 5
# stats: should correlation("corr") or regression slope ("reg") be returned? Default "corr
# removenames: Variables names that are not environmental variables
#    Defaulr: c("siteid", "spid", "x", "y", "occ", "group")

SimCorrReg <- function(species, nsim=5, stats="corr", 
                       removenames=c("siteid", "spid", "x", "y", "occ", "group"), 
                       ...) {
  if(!stats%in%c("corr", "reg")) stop("stats should be 'corr' or 'reg'")
  # Get PA data to fit a model to
  PADat <- GetPAdata(species, scale = FALSE)
  PODat <- GetPOdata(species, scale = TRUE)
  EnvNames <- names(PODat)[!names(PODat)%in%removenames]
  Sim <- SimFromData(species, sigma=NULL, newdata=NULL, nsim=nsim)
  
  GetFn <- ifelse(stats=="corr", GetCorrCoef, GetRegCoef)  
  res <- sapply(1:nsim, GetCorrCoef, PAd=PADat, POd=PODat, sims=Sim, envnames=EnvNames)
  
  res
}


# Calculate CV statistics for a species, using spatial cross-validation
#  sp: species
#  k: number of folds. Default to 5
#  size: Size of blocks. Default: 5e4
# classes: classes of models "lqpth"
# Sizes: Vector of default sizes

CalcCVStatistic <- function(sp, POdat = NULL, k=5, size=5e4, classes="lqpth", Sizes=NULL, spatial=FALSE, 
                            removenames = c("siteid", "spid", "x", "y", "occ", "group", "geometry")) {
  require(blockCV)
  
  if(!is.null(Sizes)) {
    region <- toupper(gsub("[0-9]*", "", sp))
    size <- Sizes[region]
  }
  if(is.null(data))  POdat <- GetPOdata(sp=sp, scale=TRUE, sf=TRUE)
  
  if(spatial) {
    dat.CV <- cv_spatial(x = POdat,column = "occ", 
                         k = k, # number of folds
                         size = size, # size of the blocks in metres
                         selection = "random", iteration = 50, biomod2 = FALSE,
                         progress = FALSE, report = FALSE,plot = FALSE)
    
  } else {
    # dat.CV$folds_list
    dat.CV <- list(folds_list = MakeCV(x = POdat, k = k))
  }
  
  stats <- lapply(dat.CV$folds_list, function(block, alldat, removenames) {
#    block <- dat.CV$folds_list[[1]]; alldat <- POdat
    traindat <- alldat[block[[1]],, drop=TRUE]
    testdat <- alldat[block[[2]],, drop=TRUE]
    EnvNames <- names(traindat)[!names(traindat)%in%removenames]
    #  Fit MaxEnt model
    MaxNet.mod <- FitMaxEnt(dat=traindat, EnvNames=EnvNames, RespName="occ", 
                            PA=FALSE, classes=classes)
    if(is.null(MaxNet.mod)) {
      Valid <- c(Sensitivity=NA, Specificity=NA, TSS=NA, AUC=NA)
    } else {
      testdat$test.pred <- predict(MaxNet.mod, testdat, type="link")
      Valid <- CalcFitStats(pres=testdat$occ, pred=testdat$test.pred, thresh=NULL)
    }
    Valid
  }, alldat=POdat, removenames=removenames)
  
#  stats
  Mean <- c(
    AUC = mean(unlist(lapply(stats, function(x) x["AUC"]))),
    TSS = mean(unlist(lapply(stats, function(x) x["TSS"])))
  )
  Mean
}



# non.spatial CV
# might have to balance numers of rpesences & absences too
MakeCV <- function(x = POdat, k = k) {#column = "occ", 
  Folds <- rep(1:k, times=floor(nrow(x)/k))
  if(nrow(x)%%5 != 0 ) Folds <- c(Folds, 1:(nrow(x)%%5))
  Fold.rand <- sample(Folds, replace = FALSE)
  foldsl <- lapply(1:k, function(f, rand) {
    res <- list(which(rand!=f), which(rand==f))
  }, rand=Fold.rand)
  foldsl
}


FitModels <- function(Data, classes, verbose=FALSE, link="logit", 
                      #                          removenames=c("siteid", "spid", "x", "y", "occ", "group"), otherSpBG = FALSE, 
                      valid=FALSE, pred=FALSE, prob = FALSE, 
                      savemodels = FALSE, ...) {
  if(!"PO"%in%names(Data)) stop("Data should have PO")
  if(!"PA"%in%names(Data)) stop("Data should have PA")
  if(!"EnvNames"%in%names(Data)) stop("Data should have EnvNames")
  
  require(glmnet)
  if(verbose) message("Starting ", Data$species, " now")
  
  #  Fit MaxEnt model
  MaxNet.mod <- FitMaxEnt(dat=Data$PO, EnvNames=Data$EnvNames, RespName="occ", 
                          PA=FALSE, classes=classes)
  
  if(verbose) message("Maxent for ", Data$species, " done")
  
  # Validate on PA data
  if(is.null(MaxNet.mod)) {
    PAmod <- NULL
    res <- NULL
  } else {
    Data$PA$MxPred <-  predict(MaxNet.mod, Data$PA, type="link")
    Data$PA$MxPred <-  Data$PA$MxPred - mean(Data$PA$MxPred)
    validmod <- glm(PresAbs~MxPred, data=Data$PA, family=binomial(link))
    
    # Fit MaxEnt model to PA data as PA
    PAmod <- FitMaxEnt(dat=Data$PA, EnvNames=Data$EnvNames, RespName="PresAbsNum", 
                       PA=TRUE, classes=classes)
    res <- list(coefficients = coef(validmod),
                alpha = MaxNet.mod$alpha, 
                confint = confint(validmod))
    
  }
  
  if(verbose) message("Validation for ", Data$species, " done")
  # Combine predictions  
  if(pred | valid) {
    if(!is.null(PAmod)) {
      Predicted <- data.frame(maxnet = Data$PA$MxPred,
                              valid = predict(validmod, newdata=Data$PA, type = "link"),
                              PA = predict(PAmod, newdata=Data$PA[,Data$EnvNames], type = "link")
                              
      )  
      Predicted <- apply(Predicted, 2, scale, scale=FALSE) # Mean centre
      if(valid) {
        Valid <- apply(Predicted, 2, function(pred, pres)  {
          CalcFitStats(pres=pres, pred=pred, thresh=NULL)
        }, pres=Data$PA$PresAbsNum)
      }
      
    } else {
      Predicted <- NULL
      Valid <- NA
    }
    
    # Return probabilities
    if(prob) {
      if(!is.null(PAmod)) {
        Predicted.prob <- data.frame(maxnet = predict(MaxNet.mod, Data$PA, type="logistic"),
                                     valid = predict(validmod, newdata=Data$PA, type = "response"),
                                     PA = predict(PAmod, newdata=Data$PA[,Data$EnvNames], type = "logistic")
        )
      } else {
        Predicted.prob <- NULL
      }
      
    }
  }
  # Calculate validation Statistics
  
  if(pred) {
    res$pred <- Predicted
    if(prob) res$pred.prob <- Predicted.prob
  }
  if(valid) res$valid <- Valid
  if(savemodels) {
    res$MaxEnt <- MaxNet.mod
    res$valid <- validmod
    res$PAMaxEnt <- PAmod
  }
  res
}


# Get GoF statistics, either from data or simulate data for a species
GetStatistics <- function(sp, classes="l", verbose=FALSE, sim=FALSE, outfile=NULL,
                          size = 5e4, future.seed=TRUE, link="logit", ...) {
  if(sim) {
    Data <- SimStatsFromData(sp, classes=classes, ...)
  } else {
    Data <- GetData(sp, classes=classes, ...)
  }
  if(verbose) cat("Got ", sp , " data\n")
  
  FM <- FitModels(Data, classes="l", verbose=FALSE, valid=TRUE, pred = TRUE, 
                  onlySumm = TRUE,
                  future.seed=TRUE, link="logit")  
  if(verbose) cat("Fitted ", sp , " model\n")

  if(!is.null(FM)) {
    # Cross validation
    CV <- CalcCVStatistic(sp, POdat = Data$PO, k=5, size=size, classes=classes, Sizes=NULL, spatial=FALSE)
    if(verbose) cat("Done ", sp , " CV\n")
    
    res <- c(
      corMxPA =ifelse(is.null(FM$pred)|all(is.na(FM$pred)),
                             NA, cor(FM$pred[,"maxnet"], FM$pred[,"PA"])),
             AUC = ifelse(is.null(FM$valid)| all(is.na(FM$valid)), 
                          NA, FM$valid["AUC", "maxnet"]),
             TSS = ifelse(is.null(FM$valid)| all(is.na(FM$valid)), 
                          NA, FM$valid["TSS", "maxnet"]),
             beta = FM$coefficients["MxPred"], 
             AUC.CV=CV["AUC"], 
             TSS.CV=CV["TSS"])
  } else {
    res <- rep(NA, 6)
  }
  names(res) <- c("corMxPA", "AUC", "TSS", "beta", "AUC.CV", "TSS.CV")
  if(!is.null(outfile)) cat(sp, res, file=outfile, sep=",", 
                            append = TRUE, fill=TRUE)
  res
}


# Get Gof stats from simulated data
RunSimulatedStats <- function(region, nsims=5, sizes = 5e4, classes="l", 
                              POthreshold = 20, ...) {
  size <- ifelse(length(sizes)>1 & region%in%names(sizes), sizes[region], sizes[1]) 
  
  Spp <- table(disPo(region)$spid)
  SpNames <- names(Spp)[Spp>=POthreshold]
  
  SpeciesToRun <- rep(SpNames, each=nsims)
  SimStats <- future_sapply(SpeciesToRun, GetStatistics, classes=classes, verbose=FALSE, sim=TRUE, sigma=NULL,
                            size=size, link="logit", onlySumm = TRUE, ...)
  # SimStats <- sapply(SpeciesToRun, GetStatistics, classes=classes, verbose=FALSE, sim=TRUE, sigma=NULL,
  #                    size=size, link="logit", onlySumm = TRUE, ...)
   # SimStats <- sapply(SpeciesToRun, GetStatistics, classes=classes, verbose=FALSE, sim=TRUE, sigma=NULL,
   #                    size=size, link="logit", onlySumm = TRUE)
  SimStats.df <- as.data.frame(t(SimStats))
  SimStats.df$species <- colnames(SimStats)
#  if(!is.null(outfile)) save(SimStats.df, file=paste0("SimStats", region, ".RData"))
  save(SimStats.df, file=paste0("Results/SimStats", region, ".RData"))
  SimStats.df
  
}





# GoF statistics
ExtractGof <- function(lst.sp, model="maxnet") {
  c(beta = as.numeric(lst.sp$coefficients["MxPred"]), AUC = lst.sp$valid["AUC", model],
    TSS = lst.sp$valid["TSS", model])
}


PlotGoF <- function(nm, lst, statx, staty, ...) {
  df <- lst[[nm]]
  plot(df[,statx], df[,staty], type="n", xlab="", ylab="", main=nm,...) 
  if(staty=="beta") rect(-100, 0, 100, 1, col="pink", border=NA)
  if(statx=="AUC") abline(v=0.7, lty=2)
  if(statx=="TSS") abline(v=0.0, lty=2)
  
  points(df[,statx], df[,staty]) 
}

PlotAUC <- function(nm, lst, AtRho=c(0,0), AUC=FALSE, ...) {
  Cor <- round(cor(lst[[nm]][,"int"], lst[[nm]][,"ext"], use="pairwise.complete.obs"), 2)
  plot(lst[[nm]][,"int"], lst[[nm]][,"ext"], type="n", xlab="", ylab="", main=nm, ...) 
  if(AUC) {
    rect(-1, -1, 0.5, 2, col=rgb(1,0,0,alpha=0.5), border=NA) 
    rect(-1, -1, 2, 0.5, col=rgb(1,0,0,alpha=0.5), border=NA) 
  }
  text(AtRho[1], AtRho[2], bquote(rho == .(Cor)))
  abline(0,1, lty=3)
  points(lst[[nm]][,"int"], lst[[nm]][,"ext"]) 
}

# Predictions & Correlations

GetPreds <- function(results) {
  preds <- lapply(results, function(lst) {
    cbind(lst$pred[,"PA"], lst$pred[,"valid"])
  })
}
CalcCorrs <- function(lst) {
  unlist(lapply(lst, function(x) cor(x)[1,2]))
}
PlotCorrs <- function(corrs) {
  plot(corrs$cor, jitter(as.numeric(corrs$region)), yaxt="n", xlim=c(-1,1),
       xlab="Correlation", ylab="")
  axis(2, at=1:6, labels = levels(corrs$region), las=1)
  abline(v=0, lty=3)
}



# Get numbers of presences in PO & PA
GetNPres <- function(sp, ...) {
  PA <- GetPAdata(sp, ...)
  nPresPA <- sum(PA$PresAbs==1)
  nPA <- nrow(PA)
  
  PO <- GetPOdata(sp, ...)
  nPresPO <- sum(PO$occ==1)
  nPO <- nrow(PO)
  
  c(nPA=nPresPA, totPA=nPA, nPO=nPresPO, totPO=nPO)
}
