# sp <- BBSData$species[1,]; GBIF <- gbifbirds.sf; climStack <- climStack; pseudos=Pseudos; minN =25; maxN=1e4; trainProp=0.8




FitMaxEntJAGS <- function(sp, GBIF, climStack, pseudos, minN =25, maxN=1e4, trainProp=0.8, 
                          NChains=2, nburn=1e1, niter=1e2, 
                      saveFolder=NULL, returnRes =TRUE) {
  if(is.null(saveFolder) & !returnRes) stop("This will not return or save anything")
  require(maxnet)
  gbifOccs <- GBIF[GBIF$species == unlist(sp["sciname"]),]
  if(!(nrow(gbifOccs) > minN)) {
    warning(paste0(GBIF$species, " has less than 26 observations in the GBIF data, skipping"))
    return(NULL)
  } else {
    gbifOccs <- unique(gbifOccs, .keep_all= TRUE) #Keep only unique lat/long observations
    
    #If there are more than 10,000 records that would be include in training dataset (12,500 * 80%) keep that many at random
    if(nrow(gbifOccs) > maxN/trainProp){
      gbifOccs <- gbifOccs[sample(nrow(gbifOccs), maxN/trainProp, replace = FALSE), ]
    }
    UseInTesting <- sample(nrow(gbifOccs), nrow(gbifOccs)*(1-trainProp), replace = FALSE)
    gbiftesting <- gbifOccs[UseInTesting,]
    gbiftraining <- gbifOccs[-UseInTesting,]
    
    #Run the maxent model, use MaxNet
    Data <- rbind(extract(x=climStack, y=gbiftraining),
                  extract(x=climStack, y=pseudos))[,-1]
    
    PresPseudo <- rep(c(1,0), times=c(nrow(gbiftraining), nrow(pseudos)))
    Remove <- apply(Data, 1, function(v) !any(is.na(v)))
    Data <- Data[Remove,]; PresPseudo <- PresPseudo[Remove]
    
#    myModel <- maxnet(p=PresPseudo, data=Data)
    
# start maxnet code...    
    form <- maxnet.formula(PresPseudo, Data) 
    mm <- model.matrix(form, Data)
    reg <- maxnet.default.regularization(PresPseudo, mm) * 1
    weights <- PresPseudo + (1 - PresPseudo) * 100
    lambda <- 10^(seq(4, 0, length.out = 200)) * sum(reg)/length(reg) * 
      sum(PresPseudo)/sum(weights)
    
    UseInMCMC <- c(which(PresPseudo==1), 
                   sample(which(PresPseudo!=1), size=5*sum(PresPseudo), replace = FALSE))
    
    data.jags <- list(y = PresPseudo, W = rep(1, length(PresPseudo)),
                       X = mm,
                       N = nrow(mm), NFeat = ncol(mm), 
                      MeanLambda = mean(lambda), reg = reg)
    
    useSmall <- c(sample(which(PresPseudo==0), 100),
                         sample(which(PresPseudo==1), 100))
                         
    data.jags.small <- list(y = PresPseudoSmall[useSmall], W = rep(1, length(useSmall)),
                            X = mm[useSmall,],
                            N = length(useSmall), NFeat = ncol(mm), 
                            MeanLambda = mean(lambda), reg = reg)
    
      
    MaxJags.init <- function(data) {
      p <- mean(data$y)
      lambda <- rexp(1, 1/data$MeanLambda)
      alpha <- rnorm(1, log(p/(1-p)), 0.1)
      beta <- rnorm(data$NFeat, 0, sd=(lambda*data$reg))
      res <- list(lambda =lambda, alpha = alpha, beta = beta)
      res
    }
    
    Inits <- replicate(NChains, MaxJags.init(data.jags), simplify = FALSE)
    
    library(R2jags)
    
    Params <- c("alpha", "beta")
    #Run the model 
    Model <- jags(model.file = "MaxNet.jag", 
                           data = data.jags, inits = NULL,
                           parameters.to.save = Params,
                           n.chains = NChains, n.iter = niter, n.thin = 1,
                           n.burnin = nburn)
    
    
    
    
    
    
      bb <- model$beta[, 200]
      model$betas <- bb[bb != 0]
      model$alpha <- 0
      rr <- predict.maxnet(model, data[p == 0, , drop = FALSE], 
                           type = "exponent", clamp = F)
      raw <- rr/sum(rr)
      model$entropy <- -sum(raw * log(raw))
      model$alpha <- -log(sum(rr))
      model$penalty.factor <- reg
      model$featuremins <- apply(mm, 2, min)
      model$featuremaxs <- apply(mm, 2, max)
      vv <- (sapply(data, class) != "factor")
      model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
      model$varmax <- apply(data[, vv, drop = FALSE], 2, max)
      means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
      majorities <- sapply(names(data)[!vv], function(n) which.max(table(data[p == 
                                                                                1, n, drop = FALSE])))
      names(majorities) <- names(data)[!vv]
      model$samplemeans <- unlist(c(means, majorities))
      model$levels <- lapply(data, levels)
      model
    }
    
    
    
    #Create the raw output from the spatial prediction 
    rawOut <- predict(myModel, as.data.frame(climStack), type="link")
    
    res <- list(model=myModel, predictions=rawOut, 
                presence=gbifOccs, testing=UseInTesting)
    if(!is.null(saveFolder)) {
      save(res, file=file.path(saveFolder, paste0(sp["sciname"], "Res.RData")))
    }
    if(!returnRes) res <- NULL
    res
  }
}
