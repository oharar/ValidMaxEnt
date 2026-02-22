library(maxnet)
library(disdat)
library(future.apply)
library(sdm)
library(Metrics)

source("MaxEntBayesFunctions.R")
RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")


plan(multisession, workers = 4) ## Run in parallel on local computer

# plan(sequential) ## Run in parallel on local computer

AllCoefs.lqpt <- sapply(c("AWT", "CAN", "NSW", "NZ", "SA", "SWI"), JustMaxEnt, 
                   remove=RemoveNames, classes="lqpt")

AllCoefs.l <- sapply(c("AWT", "CAN", "NSW", "NZ", "SA", "SWI"), JustMaxEnt, 
                   remove=RemoveNames, valid = TRUE, 
                   pred = TRUE, classes="l", simplify=FALSE)

# bgEnv <- disBg("CAN")
# EnvNames <- names(bgEnv)[!(names(bgEnv)%in%RemoveNames)]
# fit02 <- FitMaxEntToSp(sp="can02", dataname = "CAN", EnvNames=EnvNames, classes="l", 
#               verbose = TRUE, valid = TRUE, 
#               pred = TRUE, link = "logit") 


# CanCoefs <- JustMaxEnt(dataname="CAN", remove=RemoveNames, classes="lqpht") # runs 
# CanCoefscll <- JustMaxEnt(dataname="CAN", remove=RemoveNames, classes="lqpht", link="cloglog") # runs 
# 
# par(mfrow=c(1,2))
# plot(t(CanCoefs), main="CAN"); abline(h=1); abline(h=0)
# plot(t(CanCoefscll), main="CAN"); abline(h=1); abline(h=0)

PlotCoefs <- function(nm, lst) {
  coefs.l <- lapply(lst[[nm]], function(x) x$coeficients)
  coefs <- matrix(unlist(coefs.l), ncol = 2, byrow = TRUE)
  colnames(coefs) <- names(lst[[nm]][[2]]$coeficients)
  rownames(coefs) <- names(lst[[nm]])
  
  plot(coefs, ylim=range(c(0,1, coefs[,"Pred"])), type="n",
       xlab="", ylab="", main=nm) 
  rect(-100, 0, 100, 1, col="pink", border=NA)
  points(coefs)
  box()
}

par(mfrow=c(2,3), mar=c(3,2,3,1), oma=c(2,2,0,0))
sapply(names(AllCoefs.l), PlotCoefs, lst=AllCoefs.l)
mtext("Intercept", 1, outer=TRUE)
mtext("Slope", 2, outer=TRUE, line=0.5)




# can18 has highest slope for Canada
# can14 has lowest slope for Canada (-0.5)
# can06 looks OK
RemoveNames <- c("siteid", "spid", "x", "y", "occ", "group")
 RemoveNames <- c(RemoveNames, "ontveg")

CanPAPA <- disPa("CAN"); CanPAEnv <- disEnv("CAN")
CanPA <- merge(CanPAPA, CanPAEnv)
CanPA$ontveg <- relevel(factor(CanPA$ontveg), ref = "4")
CanPA[,grep("can", names(CanPA))] <- CanPA[,grep("can", names(CanPA))]==1

CanPO <- disPo("CAN"); CanBG <- disBg("CAN")
EnvNames <- names(CanBG)[!(names(CanBG)%in%RemoveNames)]
Dat14 <- rbind(CanPO[CanPO$spid=="can14",], CanBG)
Dat14$ontveg <- relevel(factor(Dat14$ontveg), ref = "4")


MaxNet.mod <- maxnet(p=Dat14$occ, data=Dat14[,EnvNames],
                     f=maxnet.formula(p=Dat14$occ, data=Dat14[,EnvNames],
                                      classes="l"))
beta.Mx <- MaxNet.mod$beta[,200]

form <- formula(paste0("can14 ~ ", paste(EnvNames, collapse=" + ")))
PA.mod <- glm(form, data=CanPA, family="binomial")

summary(PA.mod)

names(beta.Mx) <- gsub("categorical\\(", "", names(beta.Mx))
names(beta.Mx) <- gsub("\\):", "", names(beta.Mx))

beta.Mx <- beta.Mx[names(beta.Mx)%in%names(coef(PA.mod))]
beta.glm <- coef(PA.mod)[names(coef(PA.mod))%in%names(beta.Mx)]

plot(beta.Mx, beta.glm, type="n") #, ylim=c(-2,2))
text(beta.Mx, beta.glm , labels=names(beta.glm))

table(CanPA$can14, CanPA$ontveg)
table(Dat14$occ, Dat14$ontveg)

Pred.Mx <-  predict(MaxNet.mod, CanPA, type="link")
Pred.glm <-  predict(PA.mod, CanPA, type="link")

plot(Pred.Mx, Pred.glm)
cor(Pred.Mx, Pred.glm)
boxplot(Pred.Mx ~ CanPA$can14)

CalcFitStats <- function(pres, pred, thresh=NULL) {
  require(Metrics)
  if(is.null(thresh)) thresh <- mean(pred)
  PredOne <- pred > thresh
  
  ConfMat <- table(pres, PredOne)
  Sens <- sum(PredOne*pres)/sum(PredOne)
  Spec <- sum((1-PredOne)*(1-pres))/sum(1-PredOne)
  TSS <- Sens + Spec - 1
  AUC <- Metrics::auc(CanPA$can14, Pred.Mx)
  c(Sensitivity = Sens, Specificity = Spec, TSS=TSS, AUC=AUC)  
}

CalcFitStats(pres=CanPA$can14, pred=Pred.Mx, thresh=NULL)
CalcFitStats(pres=CanPA$can14, pred=Pred.glm, thresh=NULL)


Pred.Mx <-  predict(MaxNet.mod, CanPA, type="link")
Pred.glm <-  predict(PA.mod, CanPA, type="link")


# Try with other specie presence

sapply()
CANPO <- disPo("CAN")

POSites <- sapply(c("AWT", "CAN", "NSW", "NZ", "SA", "SWI"), function(name) {
  dat <- disPo(name)
#  table(table(dat$siteid))
  table(table(dat$y, dat$x))
})


sp="can01"; dataname="CAN"; classes="l"
  
region <- "CAN"; results <- AllCoefs.l; stat="TSS"
GetStat <- function(region, results, stat, fit=NULL) {
  gstat <- function(lst, st) lst$valid[st,]
  gcoef <- function(lst, st) lst$coefficients[st]
  
  fn <- switch(stat,
               "(Intercept)" = gcoef, 
               "Pred" = gcoef,
               "Sensitivity" = gstat, 
               "Specificity" = gstat, 
               "TSS" = gstat, 
               "AUC" = gstat, NA)
  stats <- t(list2DF(lapply(results[[region]], fn, st=stat)))
  if(ncol(stats)==3)  colnames(stats) <- colnames(results[[region]][[1]]$valid)
  if(ncol(stats)==2)  colnames(stats) <- names(results[[region]][[1]]$coefficients)
  if(fit%in%colnames(stats)) stats <- stats[,fit] # If specified, only use set column(s)
  stats
}

TSSes <- lapply(names(AllCoefs.l), GetStat, results=AllCoefs.l, stat="TSS", fit="valid")
AUCs <- lapply(names(AllCoefs.l), GetStat, results=AllCoefs.l, stat="AUC")
Slopes <- lapply(names(AllCoefs.l), GetStat, results=AllCoefs.l, stat="Pred")
names(TSSes) <- names(Slopes) <- names(AUCs) <- names(AllCoefs.l)


plot(TSSes$CAN[,"valid"], Slopes$CAN[,1])
  
TSSes[["AWT"]][,"maxnet"]
TSSes[["CAN"]][,"valid"]


AllCoefs.l <- sapply(Regions, JustMaxEnt, remove=RemoveNames, classes="l", 
                     valid = TRUE, pred=TRUE, prob=TRUE)


PlotStats <- function(region, stat1, stat2, fit, Allcoefs) {
  stat1vals <- GetStat(region, results=Allcoefs, stat=stat1, fit)
  stat2vals <- GetStat(region, results=Allcoefs, stat=stat2, fit)
  plot(stat1vals, stat2vals, xlab=stat1, ylab=stat2, main=region, xlab="", ylab="")
}

par(mfrow=c(2,3), mar=c(3,2,4,1), oma=c(2,2,0,0))
sapply(names(AllCoefs.l), PlotStats, stat1="AUC", stat2="Pred", fit="PA", Allcoefs=AllCoefs.l)
mtext("AUC", 1, outer=TRUE)
mtext("Slope", 2, outer=TRUE)

par(mfrow=c(2,3), mar=c(3,2,4,1), oma=c(2,2,0,0))
sapply(names(AllCoefs.l), PlotStats, stat1="TSS", stat2="Pred", fit="PA", Allcoefs=AllCoefs.l)
mtext("TSS", 1, outer=TRUE)
mtext("Slope", 2, outer=TRUE)





par(mfrow=c(5,4), mar=c(2,2,1,1), oma=c(2,2,0,0))
sapply(names(AllCoefs.l$CAN), function(nm, allC) {
  lst <- allC[[nm]]
  plot(lst$pred[,"PA"], lst$pred[,"valid"], xlab="", ylab="", 
       main=nm)
}, allC=AllCoefs.l$CAN)
mtext("Fitted to PA data", 1, outer=TRUE)
mtext("Fitted to PO data", 2, outer=TRUE)




pred.l <- AllCoefs.l$CAN$can01$pred[,"PA"] + AllCoefs.l$CAN$can01$coefficients["(Intercept)"]
data <- disPa("CAN")$can01
# RankPred <- rank(pred)
Resid <- data - AllCoefs.l$CAN$can01$pred.prob$PA

hh <- hist(AllCoefs.l$CAN$can01$pred.prob$PA)


Cut <- cut(pred, breaks =hh$breaks)

meanPred <- tapply(pred, list(Cut), mean)
meanResid <- tapply(Resid, list(Cut), mean)

plot(meanPred, meanResid)
plot(pred, Resid)
plot(pred, data)


