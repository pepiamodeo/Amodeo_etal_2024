#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)

#You may wish to loop over samples_list and thinning_list as done in Script S3 
nChains = 4
thin = 100
samples = 500
print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
filename_in = paste("models/models_thin_", as.character(thin),
                    "_samples_", as.character(samples),
                    "_chains_",as.character(nChains),
                    ".Rdata",sep = "")
load(file = filename_in) #models, modelnames
nm = length(modelnames)

MF = list()
MFCV = list()
#MFCV.plot = list()
WAIC = list()

for(model in 1:nm){
  m<-get(modelnames[model])
  print(paste0("model = ",as.character(model)))
  preds = computePredictedValues(m)
  MF[[model]] = evaluateModelFit(hM=m, predY=preds)
  ####### ACTIVAR #####
  # cross val sample level (transect)
  ##partition = createPartition(m, nfolds = 5) # five fold crossvalidation 80% fit - 20% predict
  ##preds = computePredictedValues(m,partition=partition, nParallel = nChains)
  ##MFCV[[model]] = evaluateModelFit(hM=m, predY=preds)
  #####################
  # cross val plot level
  #partition = createPartition(m, nfolds = 5, column="plot") # five fold crossvalidation 80% fit - 20% predict
  #preds.plot = computePredictedValues(m,partition=partition) # fold por plots
  #MFCV.plot[[model]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[model]] = computeWAIC(m)
}

filename_out = paste("models/MF_thin_", as.character(thin),
                                            "_samples_", as.character(samples),
                                            "_chains_",as.character(nChains),
                                            ".Rdata",sep = "")
save(MF,MFCV,WAIC,modelnames,file = filename_out)

