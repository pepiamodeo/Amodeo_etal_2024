
library(plyr)
library(Hmsc)

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

# Definicion de modelo PA
# defino dos, uno con efectos fijos aditivos entre las covariates X
# otro con interacciones. para simplificar el modelo pongo bioclimate como random

# community data matrix Y,
# the environmental data frame XData,
# and the phylogenetic tree PhyloTree.

# acomodo nombre para que coincidan con este script
PData<-P
XData<-X
SData<-S
TrData<-Tr

# exploro YData
dim(Y) # 299 sampling units and 26 species

# Descarto especies raras

Y[,colSums(Y)==0] # especies con cero

colSums(Y>0) # veo en cuantas samples aparecio cada especie

spSEL = colSums(Y >0)>=10
sum(spSEL) # cantidad de especies seleccionadas

sum(Y[,spSEL])/sum(Y) # proporcion de la abundancia total que abarco

Ytot <- Y
Y <- Y[, spSEL] # aplico el filtro de seleccion de especies

#saco reporte general de la matriz total Ytot
# Genero YData para modelo PA y AvCOP

Y.PA = 1*(Ytot>0) # transformo en datos de Presencia Ausencia

# abundance COP, log transformed p 193/370 Ovaskainen Abrego 2020
Y.Ab<-Ytot
Y.Ab[Y.Ab==0] = NA
Y.Ab = log(Y.Ab)
Y.Ab<-as.matrix(Y.Ab)

# We next explore the distribution of species richness S as the row sums, and the distribution species prevalence as column means

Riq = rowSums(Y.PA) # species richness
Prev = colMeans(Y.PA) # prevalence
Arel = colSums(Ytot)/sum(Ytot) # abundance relative
Ntot = rowSums(Ytot) # total abundance

summary(Riq) # riqueza
summary(Prev) # prevalencia
summary(Arel) # abundancia relativa por sp
summary(Ntot) # abundancia absoluta por transecta

hist(Riq, xlab = "Species richness (S)")
hist(Prev, xlab = "Species prevalence (P)",xlim=c(0,1))
hist(Arel, xlab = "Relative Abundance (A)",xlim=c(0,1))

sort(Arel,decreasing = T)
sort(Prev,decreasing = T)

# exporto tabla con prevalencias y abundancias relativas totales de todo el dataset

TAB.SPP<-merge(data.frame(Sp=names(sort(Prev,decreasing = T)),P=sort(Prev,decreasing = T)),
      data.frame(Sp=names(sort(Arel,decreasing = T)),A=sort(Arel,decreasing = T)))
TAB.SPP<-TAB.SPP[order(TAB.SPP$P,decreasing = T),]
write.csv(TAB.SPP,"./panels/tab_sp.csv")

# REPITO EL REPORTE con las especies seleccionadas

# Genero YData para modelo PA y AvCOP

Y.PA = 1*(Y>0) # transformo en datos de Presencia Ausencia

# abundance COP, log transformed p 193/370 Ovaskainen Abrego 2020
Y.Ab<-Y
Y.Ab[Y.Ab==0] = NA
Y.Ab = log(Y.Ab)
Y.Ab<-as.matrix(Y.Ab)

# We next explore the distribution of species richness S as the row sums, and the distribution species prevalence as column means

Riq = rowSums(Y.PA) # species richness
Prev = colMeans(Y.PA) # prevalence
Arel = colSums(Y)/sum(Y) # abundance relative
Ntot = rowSums(Y) # total abundance

range(Riq)
range(Prev)
range(Arel)
range(Ntot)

hist(Riq, xlab = "Species richness (S)")
hist(Prev, xlab = "Species prevalence (P)",xlim=c(0,1))
hist(Arel, xlab = "Relative Abundance (A)",xlim=c(0,1))

sort(Arel,decreasing = T)
sort(Prev,decreasing = T)

# For species richness, the y-axis (Frequency) corresponds to the number of sampling units,
# and the variable in the x-axis (S) is the number of species found from each sampling unit.
# For species prevalence, the y-axis (Frequency) corresponds to the number of species,
# and the variable in the x-axis (P) is the fraction of sampling units in which the species is present.
# The prevalence of the species range from 27% to 99%, reflecting the fact that we have included only the most common species in these data.
# Species richness varies greatly among the sampling units, ranging from 0 to 47.

# Let us then look at the environmental and spatial data
head(XData)

# the following also show the types (factor, num) of variables:
str(XData)

# reordeno variables categoricas
XData$Bioclimate<-factor(XData$Bioclimate,levels=c("Supra","Meso","Termo"))
XData$Lithology<-factor(XData$Lithology,levels=c("Limestone","Alluvial","Marl","Metamorphic"))
XData$Type<-factor(XData$Type,levels=c("OS","TL","SB","NV"))
XData$Ext.SB<-as.factor(XData$Ext.SB)
XData$Ext.NV<-as.factor(XData$Ext.NV)
XData$Legacy <-as.factor(XData$Legacy)
str(XData)

# and a full summary
summary(XData)

# Let us then look at the data on species traits and phylogenetic relationships.

PData<-PData[spSEL,] # aplico seleccion de especies a TrData
TrData<-TrData[spSEL,] # aplico seleccion de especies a TrData

# Chequeo NA en Tr para las 12 sp seleccionadas
is.na(TrData) # traer spSEL de script S2
apply(X=TrData,MARGIN=1,FUN=anyNA)

# kjg tiene 5 especies con NA

TrData<-TrData[,-12] # saco kjg que tiene NA
TrData$Family <- PData$Family  # paso info de familia a TrData
levels(TrData$Sexual_System)[2]<-"monoecious"

TrData<-droplevels(TrData)

str(TrData)

# We are now ready to define the HMSC model. While the book defines three models (m.FULL, m.ENV, and m.SPACE),
# here we define the full model only. The remaining models are special cases of the FULL model and thus
# the script can be modified easily to fit those alternative models.
# Thus, the model (to be called here simply as m) includes both the environmental covariates as well as the spatial random effect of the route.

# To define a spatial random effect at the level of the route, we need to include the route id:s in the studyDesign

studyDesign = data.frame(plot = SData$PlotID, 
                         transect= SData$TransectID)

# We next define the random level object. If we would define an unstructured random effect, we would use the
# units = ... argument as we did in the fungal example of Section 7.9. As we wish to define a spatial
# random effect, we use instead the sData argument.

rL.plot = HmscRandomLevel(units=levels(studyDesign$plot))
rL.transect = HmscRandomLevel(units=levels(studyDesign$transect)) # version anterior con rl en transecta

# spatial random effects

SData_xy<-droplevels(SData[!duplicated(SData$PlotID),])

xycoords<-SData_xy[3:4]
rownames(xycoords)<-SData_xy$PlotID

rL.xy = HmscRandomLevel(sData = xycoords)

# Note that the row names of xy correspond to the units of the studyDesign. This is necessary to make
# Hmsc understand how the units of the random effect (rows of xy) correspond to the sampling units
# (rows of studyDesign, Y and XData). While in this example the dimensions of xy and studyDesign match,
# this is not necessarily the case. For example, if we would have included multiple surveys to the same sites,
# each survey would be a row of studyDesign, but sData (or units in case of unstructured random effect)
# should have each site only once.

# As environmental covariates, we include the habitat type as the categorical variable, and the second order response
# to climatic conditions. We use the poly function to defined the second order response. We further use in that
# raw = TRUE to make the model formulation consistent with the way in which predictions over environmental
# gradients are made in Hmsc.

XFormula3 = ~ Altitude + Lithology + Type + Age + Ext.SB+ Ext.NV + Legacy
#XFormula2 = ~ Altitude + Lithology + Type + Age + Ext.SB+ Ext.NV
#XFormula1 = ~ Bioclimate + Lithology + Type + Age + Ext.SB+ Ext.NV

names(TrData)

#TrFormula = ~ GrowthForm + LeafPhenology + LeafShape + RootDepth + SLA + Shade + Sexual_System + dispcat + leng + sdm + seeds + pcw + Family
TrFormula = ~ GrowthForm + RootDepth + SLA + Sexual_System + dispcat + leng + sdm + pcw
#TrFormula3 = ~ GrowthForm + LeafPhenology + LeafShape + RootDepth + SLA + Shade + Sexual_System + dispcat + leng + sdm + seeds + pcw + Family


# We are now ready to define the model. Note that in ranLevels=list(Route=rL), "Route" refers to a column name
# of the studyDesign

models = list()

models[[1]] = Hmsc(Y=Y.PA, XData = XData, XFormula=XFormula3,
           TrData = TrData,TrFormula = TrFormula,
           distr="probit", studyDesign=studyDesign,
           ranLevels=list(plot=rL.xy))

# en el libro plot esta entre comillas ("plot"=rL.plot)

models[[2]] = Hmsc(Y=Y.Ab, XData = XData, XFormula=XFormula3,
                   TrData = TrData,TrFormula = TrFormula,
                   distr="normal", studyDesign=studyDesign,
                   ranLevels=list(plot=rL.xy))

#models[[3]] = Hmsc(Y=Y.PA, XData = XData, XFormula=XFormula2,
#                   TrData = TrData,TrFormula = TrFormula,
#                    distr="probit", studyDesign=studyDesign,
#                    ranLevels=list(plot=rL.xy))
# 
# models[[4]] = Hmsc(Y=Y.Ab, XData = XData, XFormula=XFormula2,
#                    TrData = TrData,TrFormula = TrFormula,
#                    distr="normal", studyDesign=studyDesign,
#                    ranLevels=list(plot=rL.xy))

# As always, it is a good idea to explore the model object

m<-models[[1]] # cambiar manualmente de modelo
m
# This should give "Hmsc object with 137 sampling units, 50 species, 7 covariates, 4 traits and 1 random levels"

head(m$X)

# The factor of habitat type has been expanded to dummy variables, and the model includes first and second orders of climate

head(m$XScaled)

# In the computations, scaled versions of climate will be used

head(m$Tr)

# The factor of migratory strategy has been expanded to dummy variables, and the model includes the linear effect of log(mass)

head(m$TrScaled)

#  datos campos

table(XData$Bioclimate,XData$Lithology)
tab.campos<-cbind(studyDesign,XData)
ddply(tab.campos,.(Bioclimate,Lithology),summarise,length(unique(plot)))

write.csv(table(XData$Bioclimate,XData$Lithology),
          "./panels/tab_campos.csv")

# grabo modelos no ajustados
modelnames<-c("PA_model","AbCOP_model")
save(models,modelnames,file = "models/unfitted_models") #models, modelnames
