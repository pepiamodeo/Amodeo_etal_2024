#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)
library(colorspace)
library(vioplot)

#include in samples_list and thin_list only those models that you have actually fitted!
samples_list = c(500)
thin_list = c(100)
nst = length(thin_list)
nChains = 4

ma = NULL
na = NULL
for (Lst in 1:nst) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
  load(filename)
  nm = length(modelnames)
  #nm = 2 # selecciono modelos sin interacción
  for(j in 1:nm){
    m<-get(modelnames[j])
    mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0(as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      if(j==1){
        na = c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na = c(na,"")
      }
    }
  }
}

# resumen de los scale reduction factor
summary(ma[,1])
summary(ma[,2])
summary(ma[,3])
summary(ma[,4])

table(ma[,1]>1.10) # cuantos por encima de 1.1?
ma[,1][ma[,1]>1.10] # cuales son?


pdf(file=paste("./panels/MCMC_convergence.pdf"))
par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
dev.off()

