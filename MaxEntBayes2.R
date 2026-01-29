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

NZSpecies <- unique(disPo("NZ")$spid)
Results <- as.data.frame(t(sapply(NZSpecies, GetResults, dataname="NZ")))


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

NZ.int <- GetStat(dataname="NZ", stat="alpha")
AWT.int <- GetStat(dataname="AWT", stat="alpha")
CAN.int <- GetStat(dataname="CAN", stat="alpha")
NSW.int <- GetStat(dataname="NSW", stat="alpha")
SA.int <- GetStat(dataname="SA", stat="alpha")
SWI.int <- GetStat(dataname="SWI", stat="alpha")

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
