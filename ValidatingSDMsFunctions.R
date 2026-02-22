
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
  Sens <- sum(PredOne*pres)/sum(PredOne)
  Spec <- sum((1-PredOne)*(1-pres))/sum(1-PredOne)
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
  
  # Validate on PA data
  # Get PA data and centre
  PAdata <- GetPAdata(sp=sp, scale=FALSE)
  PAdata$p.p <- as.numeric(as.character(PAdata$PresAbs)) # convert PA to numeric from logical
  
  # better to pass these into GetPAdata()?
  if(!is.null(attr(POdata, "means"))) {
    PAdata[,names(attr(POdata, "envmeans"))] <- sweep(PAdata[,names(attr(POdata, "envmeans"))], 
                                                      2, attr(POdata, "envmeans"), "-")
    PAdata[,names(attr(POdata, "envsds"))] <- sweep(PAdata[,names(attr(POdata, "envsds"))], 
                                                    2, attr(POdata, "envsds"), "/")
  }
  
  # Validate MaxEnt model
  
  PAdata$MxPred <-  predict(MaxNet.mod, PAdata, type="link")
  PAdata$MxPred <-  PAdata$MxPred - mean(PAdata$MxPred)
  validmod <- glm(PresAbs~MxPred, data=PAdata, family=binomial(link))
  
  # Fit MaxEnt model to PA data as PA
  PAmod <- FitMaxEnt(dat=PAdata, EnvNames=EnvNames, RespName="p.p", 
                     PA=TRUE, classes=classes)
  
  if(verbose) message("Validation for ", sp, " done")
  # Combine predictions  
  if(pred | valid) {
    if(!is.null(PAmod)) {
      Predicted <- data.frame(maxnet = PAdata$MxPred,
                              valid = predict(validmod, newdata=PAdata, type = "link"),
                              PA = predict(PAmod, newdata=PAdata[,EnvNames], type = "link")
      )  
    } else {
      Predicted <- NULL
    }
    
    # Return probabilities
    if(prob) {
      if(!is.null(PAmod)) {
        Predicted.prob <- data.frame(maxnet = predict(MaxNet.mod, PAdata, type="logistic"),
                                     valid = predict(validmod, newdata=PAdata, type = "response"),
                                     PA = predict(PAmod, newdata=PAdata[,EnvNames], type = "logistic")
        )
      } else {
        Predicted.prob <- NULL
      }
      
    }
    Predicted <- apply(Predicted, 2, scale, scale=FALSE) # Mean centre
  }
  # Calculate validation Statistics
  if(valid) {
    Valid <- apply(Predicted, 2, function(pred, pres)  {
      CalcFitStats(pres=pres, pred=pred, thresh=NULL)
    }, pres=PAdata$p.p)
  }
  res <- list(coefficients = coef(validmod),
              alpha = MaxNet.mod$alpha, 
              confint = confint(validmod))
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


# Just fit & validate maxEnt
JustMaxEnt <- function(region, remove, classes="l", verbose=FALSE, 
                       link="logit", ...) {
  require(disdat)
  require(maxnet)
  require(future.apply)
  bgEnv <- disBg(region)
  EnvNames <- names(bgEnv)[!(names(bgEnv)%in%remove)]
  SpNames <- unique(disPo(region)$spid)
  
  #  sp <- SpNames[38]
  #  if(region=="SWI") classes <-"l"
  Coefs <- future_sapply(SpNames, FitMaxEntToSp, 
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
  rect(-100, 0, 100, 1, col="pink", border=NA)
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

GetPAdata <- function(sp, removenames = c("siteid", "spid", "x", "y", "occ", "group"), 
                      scale=TRUE) {
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
  
  if(scale) {
    EnvNamesNoFactor <- names(PAEnv)[!names(PAEnv)%in%removenames & !sapply(PAEnv, is.factor)]
    means <- apply(PAEnv[,EnvNamesNoFactor], 2, mean)
    sds <- apply(PAEnv[,EnvNamesNoFactor], 2, sd)
    PAEnv[,EnvNamesNoFactor] <- scale(PAEnv[,EnvNamesNoFactor])
    attr(PAEnv, "envmeans") <- means
    attr(PAEnv, "envsds") <- sds
    
  }
  
  PAdata <- merge(PA, PAEnv, by="siteid")
  PAdata
}


# Get PO data
# Returns data frame with "occ, a 0/1 integeter: 0=background, 1=presence, 
#  and columns of the environmental coavriates, converted to factors when appropriate.
GetPOdata <- function(sp, removenames = c("siteid", "spid", "x", "y", "occ", "group"), 
                      otherSpBG=FALSE, scale=TRUE) {
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
  Pres <- disPo(region)[disPo(region)$spid==sp,c("occ", EnvNames)]
  
  dat <- rbind(Pres, bgEnv[,c("occ", EnvNames)])
  dat <- ConvertFactors(dat, region)
  
  if(scale) {
    EnvNamesNoFactor <- names(dat)[!names(dat)%in%removenames & !sapply(dat, is.factor)]
    means <- apply(dat[,EnvNamesNoFactor], 2, mean)
    sds <- apply(dat[,EnvNamesNoFactor], 2, sd)
    dat[,EnvNamesNoFactor] <- scale(dat[,EnvNamesNoFactor])
    attr(dat, "envmeans") <- means
    attr(dat, "envsds") <- sds
  }
  dat
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
    SimPO <- sample.int(n=length(Wt), size=N, prob=Wt/sum(Wt))
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
  POData <- GetPOdata(sp=species, scale=FALSE)
  N <- sum(POData$occ)
  if(is.null(newdata) | PA) {
    PAData <- GetPAdata(sp=species, scale=FALSE)
  }
  if(PA) {
    Data <- PAData
  } else {
    Data <- POData
  }
  EnvNames <- names(Data)[!(names(Data)%in%removenames)]
  
  # fit model with linear features to data. If PO, uses (almost) infinite weights.
  YY <- ifelse(PA, "PresAbs ~", "occ ~")
  f <- formula(paste(YY, paste(EnvNames, collapse = " + ")))
  
  if(!PA) {
    Wt <- 1+ (Data$occ * 99)
  } else {
    Wt <- rep(1, nrow(Data))
  }
  truemod <- glm(f, data=Data, family=binomial("logit"), weights=Wt)
  
  # Simulate data
  # if no new data, simulate on orib´ginal PA and PO data
  if(is.null(newdata)) {
    PAData$eta <- predict(truemod, newdata = PAData, type = "link")
    POData$eta <- predict(truemod, newdata = POData, type = "link")
    
    SimPODat <- replicate(nsim, SimDataFromWeights(lnWt=POData$eta, PA=FALSE, N=N, sigma=sigma[1]))
    SimPADat <- replicate(nsim, SimDataFromWeights(lnWt=PAData$eta, PA=TRUE, sigma=sigma[2]))
  } else {
    newdata$eta <- predict(truemod, newdata = newdata, type = "link")
    
    SimPODat <- replicate(nsim, SimDataFromWeights(lnWt=newdata$eta, PA=FALSE, N=N, sigma=sigma[1]))
    SimPADat <- replicate(nsim, SimDataFromWeights(lnWt=newdata$eta, PA=TRUE, sigma=sigma[2]))
    
  }
  
  # Simulate data
  
  res <- list(PA=PA, simPA=SimPADat, simPO=SimPODat)
  if(verbose) {
    res$truemod=truemod
    res$data=newdata
  }
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
  
  PADat <- GetPAdata(species, scale = FALSE)
  PODat <- GetPOdata(species, scale = TRUE)
  EnvNames <- names(PODat)[!names(PODat)%in%removenames]
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
    mod <- tryCatch(maxnet(p=dat[,RespName], data=dat[,EnvNames], wt=wt, 
                           f=maxnet.formula(p=dat[,RespName], data=dat[,EnvNames],
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

