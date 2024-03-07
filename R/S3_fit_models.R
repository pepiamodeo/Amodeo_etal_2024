
library(Hmsc)

load(file = "models/unfitted_models") #models, modelnames

# We will store 100 posterior samples for each of two chains
# We note that for more "final" results, one might wish to have e.g. 1000 samples for each of four chains
# After fitting all models, we save the model object to a file
# Loading the fitted model then serves as the starting point for exploring the results
# The script runs over a loop where thin is first 1, then 10, then 100, and so on
# Thin is the thinning parameter of the MCMC chain.
# The transient (also called burn-in) is set to 50*thin
# When thin = 1, there will be 50 burn-in and 100 actual iterations. All actual iterations are stored.
# When thin = 10, there will be 500 burn-in and 1000 actual iterations. The actual iterations are thinned by 10, so 100 are stored.
# When thin = 100, there will be 5000 burn-in and 10000 actual iterations. The actual iterations are thinned by 100, so 100 are stored.
# A long MCMC chain is needed to achieve convergence
# Thinning is applied to avoid storing model objects of very large size
# Even if in the end thin = 1000 is required to achieve converge, We recommend to run the loop thin = 1, 10, 100, 1000
# This is for several reasons.
# First of all, it will not be known beforehand how much thinning is needed to achieve satisfactory convergence
# Second, thin = 1 will run very fast, whereas thin = 1000 will take very long (1000 times longer)
# After thin = 1 is completed, it is already possible to develop all the remaining scripts that explore the fitted model
# When exploring the fitted model, often one realizes changes that need to be made, even if the fitting has not converged
# Third, running the model fitting for thin = 1, 10, 100, 1000 does not take much longer than running it just for thin = 1000 (it takes ca. 12% longer)
# Thus, in summary, running the model fitting for thin = 1, 10, 100, 1000 typically saves a lot of time,
# as it allows one to proceed fast in writing (and revising) all the scripts that are needed from defining the model to producing the result tables and figures
# The idea is not to run the entire loop in one go, as that would take a lot of time. Just run thin = 1, and then move to develop the next scripts.
# You may then leave the remaining part of the loop (e.g. thin = 10, 100, 1000) to run e.g. overnight

# sampleo de prueba

nChains = 2
nParallel = 2
samples = 10
thin = 1
transient = 5

nm = length(models)
for (model in 1:nm){
  m <- sampleMcmc(models[[model]], thin = thin, samples = samples, 
                             transient = transient,
               nChains = nChains,
               nParallel = nParallel)
  assign(x=modelnames[model],
         value=m,
         envir = .GlobalEnv)
}

filename = paste("models/models_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")

modelnames # copiar y pegar los nombres de los objetos

#get(modelnames[1])

save(PA_model,AbCOP_model,
     modelnames,
     file=filename)

# MCMC convergence can be difficult to achieve especially in those models that are not based on normal distribution.
# For this reason, in the script above we initialize model with
# initPar="fixed effects", with which option the MCMC chains are not started from locations randomized from the prior
# but from a maximum likelihood solution to the fixed-effects part of the model

# sampleo definitivo
# definir parametros a ajustar
samples_list = 500
thin_list = c(1,10,100,1000)
nChains = 4
nParallel = 4

time.init<-Sys.time()
for(Lst in 1:length(thin_list)){
  thin = thin_list[Lst]
  samples = samples_list
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  
  for (model in 1:nm) {
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    if(model==1){initPar="fixed effects"}else{initPar=NULL}
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains,
                   nParallel = nParallel,
                   initPar = initPar) 
    assign(x=modelnames[model],
           value=m,
           envir = .GlobalEnv)
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(PA_model,AbCOP_model,
       modelnames,
       file=filename)
}
time.final<-Sys.time()
duration<-time.final-time.init
#6.69hs duracion
