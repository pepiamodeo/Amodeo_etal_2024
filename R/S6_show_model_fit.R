#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)

# NO LIMPIAR ENVIRONMENT ENTRE S5 Y S6

#This script shows model fit for probit and linear models.
#For Poisson models, you may plot e.g. pseudoR2, see e.g. the book for examples

thin = 100
samples = 500
nChains = 4

filename = paste("models/MF_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(MF)
fileout = paste("./panels/model_fit.pdf")
pdf(file = fileout)

for(j in 1:nm){
  cMF = MF[[j]]
  cWAIC = WAIC[[j]]
  m<-get(modelnames[j])
  spNames<-m$spNames
  if(!is.null(cMF$TjurR2)){
    plot(as.factor(spNames),cMF$TjurR2,
         main=paste0("explanatory power",modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": Tjur R2"),
         sub=paste("Mean: ",round(mean(cMF$TjurR2),digits=2),
                   "WAIC: ",round(cWAIC,digits=2))) # PEPI
    }
  if(!is.null(cMF$R2)){
    plot(as.factor(spNames),cMF$R2,
         main=paste0("explanatory power",modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": R2"),
         sub= paste("Mean: ",round(mean(cMF$R2),digits=2),
                    "WAIC: ",round(cWAIC,digits=2))) # PEPI
    }
  if(!is.null(cMF$AUC)){
    plot(as.factor(spNames),cMF$AUC,
         main=paste0( "explanatory power",modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": AUC"),
         sub= paste("Mean: ",round(mean(cMF$AUC),digits=2),
                    "WAIC: ",round(cWAIC,digits=2))) # PEPI
  }
  
}
dev.off()

# summary

modelnames[1]
summary(MF[[1]]$TjurR2) # explanatory power
summary(MFCV[[1]]$TjurR2) # predictive power

modelnames[2]
summary(MF[[2]]$R2) # explanatory power
summary(MFCV[[2]]$R2) # predictive power

# loop original
 for(j in 1:nm){
  cMF = MF[[j]]
  cMFCV = MFCV[[j]]
  if(!is.null(cMF$TjurR2)){
    plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": Tjur R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$R2)){
    plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$AUC)){
    plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": AUC"))
    abline(0,1)
    abline(v=0.5)
    abline(h=0.5)
  }
} ### loop original del curso power vs predictive
dev.off()

