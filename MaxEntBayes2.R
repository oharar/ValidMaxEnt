#library(BayesMaxEnt)
#library(maxnet)
library(disdat)
library(coda)

#library(INLA)
#library(parallel)


source("MaxEntBayesFunctions.R")

CAN <- disData("CAN")
EnvPO <- c("alt", "asp2", "ontprec", "ontprec4", "ontprecsd", "ontslp", 
           "onttemp", "onttempsd", "onttmin4", "ontveg", "watdist")
Species <- "can02"
CanSpecies <- unique(disPo("CAN")$spid)


# sp <- "can01"; dataname <- "CAN"
Results <- as.data.frame(t(sapply(CanSpecies, GetResults, dataname="CAN")))

# Get wide prior models
ResultsWide <- as.data.frame(t(sapply(CanSpecies, GetResults, dataname="CAN",
                                      fileprefix = "Results/MaxEndBigLambda",
                        simplify = TRUE)))
Results <- ResultsWide

# NZSpecies <- unique(disPo("NZ")$spid)
# Results <- t(sapply(NZSpecies, GetResults, dataname="NZ"))


Intercept <- Results[,grep("Interc", names(Results))]
names(Intercept) <- gsub("(Intercept).", "", names(Intercept), fixed = TRUE)
Intercept <- cbind(Intercept, 
                    CalcLims(df=Intercept, modname="MaxEnt"),
                    CalcLims(df=Intercept, modname="BayesMaxEnt"))


Slope <- Results[,!grepl("Interc", names(Results))]
names(Slope) <- gsub("\\..*\\.", ".", names(Slope))
Slope <- cbind(Slope, 
                   CalcLims(df=Slope, modname="MaxEnt"),
                   CalcLims(df=Slope, modname="BayesMaxEnt"))

PlotPost <- function(Par, adj=0.3, ...) {
  plot(Par$MaxEnt.mean,1:nrow(Par), 
       xlim=range(c(1,Par[,grep("Lim", names(Par))])), pch=16, ...)
  points(Par$BayesMaxEnt.mean, adj+(1:nrow(Par)), 
         col=2, pch=16)
  segments(Par$MaxEnt.LLim, 1:nrow(Par), 
           Par$MaxEnt.ULim, 1:nrow(Par))
  segments(Par$BayesMaxEnt.LLim, adj + (1:nrow(Par)), 
           Par$BayesMaxEnt.ULim, adj + (1:nrow(Par)), col=2)
}


par(mfrow=c(1,2), mar=c(4.1,1,1,1), oma=c(0,4,0,0))
PlotPost(Intercept, adj=0.3, xlab="Intercept", ylab="", yaxt="n")
axis(2, at=1:nrow(Intercept), labels = rownames(Intercept), las=1)  
abline(v=0)
legend(-9, 13, c("MaxEnt", "Bayes"), fill=1:2, cex=0.8)
PlotPost(Slope, adj=0.3, xlab="Slope", ylab="", yaxt="n")
#axis(2, at=1:nrow(Slope), labels = rownames(Slope), las=1)  
abline(v=0); abline(v=-1, lty=3)
mtext("Species", 2, outer=TRUE, line=2.8)


# Load MCMC chains, extract posterior summaries for alpha
# Load INLA models, extract posterior summaries for alpha
MCMC1 <- load("Results/MaxEntResMCMCNZ_l_nz01.RData")
inla1 <- load("Results/MaxEntResINLANZ_l_nz01.RData")

summary(as.mcmc.list(Fit$mcmc))
(Valid$BayesMaxEnt)

NZSpecies <-  unique(disPo("NZ")$spid)


# Fix this if we want the slope?
GetStat=function(dataname, stat) {
  Species <- unique(disPo(dataname)$spid)
  if(stat=="alpha") v.stat <- "LP" else v.stat=stat
  
  alpha <- t(sapply(Species, function(sp, dname, stat, v.stat) {
    MCMC1 <- load(paste0("Results/MaxEntResMCMC", dname, "_l_", sp, ".RData"))
    inla1 <- load(paste0("Results/MaxEntResINLA", dname, "_l_", sp, ".RData"))
    fitQuants <- summary(as.mcmc.list(Fit$mcmc))$statistics[stat, c("Mean", "SD")]
    validQuants <- unlist(Valid$BayesMaxEnt[v.stat, c("mean", "sd")])
    c(fitQuants, validQuants)
  }, dname=dataname, stat=stat, v.stat=v.stat))
  colnames(alpha) <- c("Fit.mean", "Fit.sd", "Valid.mean", "Valid.sd")
  alpha
}

IntStats <- list(
  NZ = GetStat(dataname="NZ", stat="alpha"),
  AWT = GetStat(dataname="AWT", stat="alpha"),
  CAN = GetStat(dataname="CAN", stat="alpha"),
  NSW = GetStat(dataname="NSW", stat="alpha"),
  SA = GetStat(dataname="SA", stat="alpha"),
  SWI = GetStat(dataname="SWI", stat="alpha")
)

par(mfrow=c(2,3), mar=c(2,2,4,1), oma=c(2.5,2.5,0,0))
sapply(names(IntStats), function(name, intlst) {
  int <- intlst[[name]]
  At.x <- min(int[,"Fit.mean"]) + 0.7*(diff(range(int[,"Fit.mean"])))
  At.y <- min(int[,"Valid.mean"]) + 0.8*(diff(range(int[,"Valid.mean"])))
  plot(int[,"Fit.mean"], int[,"Valid.mean"], main=name, 
       xlab="", ylab="", col="grey60")
  corrr <- round(cor(int[,"Fit.mean"], int[,"Valid.mean"]),2)
  Expr <- bquote(rho== .(corrr))
  text(At.x, At.y, Expr, cex=1.3)
  cat(At.x, " " , At.y, "\n")
}, intlst = IntStats)
mtext("Intercept for PO data", 1, line=1, outer=TRUE)
mtext("Intercept for PA data", 2, line=1, outer=TRUE)



cor(NZ.int[,"Fit.mean"], NZ.int[,"Valid.mean"])
plot(NZ.int[,"Fit.mean"], NZ.int[,"Valid.mean"])

cor(AWT.int[,"Fit.mean"], AWT.int[,"Valid.mean"])
plot(AWT.int[,"Fit.mean"], AWT.int[,"Valid.mean"])

cor(CAN.int[,"Fit.mean"], CAN.int[,"Valid.mean"])
plot(CAN.int[,"Fit.mean"], CAN.int[,"Valid.mean"])

cor(NSW.int[,"Fit.mean"], NSW.int[,"Valid.mean"])
plot(NSW.int[,"Fit.mean"], NSW.int[,"Valid.mean"])

cor(SA.int[,"Fit.mean"], SA.int[,"Valid.mean"])
plot(SA.int[,"Fit.mean"], SA.int[,"Valid.mean"])

cor(SWI.int[,"Fit.mean"], SWI.int[,"Valid.mean"])
plot(SWI.int[,"Fit.mean"], SWI.int[,"Valid.mean"])




