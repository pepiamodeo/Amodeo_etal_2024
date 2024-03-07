
# libraries
library(Hmsc)
library(colorspace)
library(corrplot)
#library(writexl)

library(reshape2)
library(plyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

library(DescTools)
theme_set(theme_bw())

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} #Extract Legend function

# load files

nChains = 4
samples = 500
thin = 100

filename = paste("models/models_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(modelnames)
#nm = 2
#postBeta = getPostEstimate(m, parName="Beta")


# Variance partitioning plot model PA ##############

  m = get(modelnames[1])
  
  # grouping variables (adaptado a cada modelo)
  var.group<-attributes(m$X)$assign+1
  m$covNames
  var.group.names<-c("Intercept","Altitude", "Lithology","Type","Age","Ext.SB","Ext.NV","Legacy")
  
  # compute model with grouping variables
  VP = computeVariancePartitioning(m,group= var.group,groupnames = var.group.names) # agrupando variables
  vals = VP$vals
  # compute explanatory powers
  preds = computePredictedValues(m)
  MF = evaluateModelFit(hM=m, predY=preds)
  
  # construyo dataframe para ggplot
  df<-melt(vals)
  names(df)<-c("Variable","Species","VarP")
  df<-droplevels(df[df$Variable!="Intercept",]) # saco el intercept
  
  #eleccion del colorpalette
  #display.brewer.all(colorblindFriendly = TRUE)
  #display.brewer.pal(n=6,name="RdBu")
  #display.brewer.pal(n=2,name="Greys")
  customcolors<-c(brewer.pal(n=6,name="RdBu"),brewer.pal(n=3,name="Greys")[1:2]) # custom colors (8) a partir de brewer colorblindfriendly
  
  # Etiquetas de las variables
  levels(df$Variable) <- c("Altitude","Lithology",
                           "Transect Type","Time since Abandonment",
                           "External TE","External NV","Legacy",
                           "Random Spatial: Plot") 
  
  var.vals<-ddply(df,.(Variable),summarise,val=round(mean(VarP),digits=2)) # calculo R2 promedio
  var.labels<-apply(var.vals, 1 , paste , collapse = " (") # pego el nombre de la variable con su R2
  var.labels<-paste(var.labels,")",sep="") # agrego parentesis que falta

  # Etiquetas de las especies
  
  fitSp<-data.frame(Species=levels(df$Species),
                    fit=round(MF$TjurR2,digits=2)) # R2 para cada espeie
  Sp.labels<-apply(fitSp, 1 , paste , collapse = " (") # pego el nombre de la especie con su R2
  Sp.labels<-paste(Sp.labels,")",sep="") # agrego parentesis que falta
  
  # gráfico
  p1<-ggplot(data=df,aes(x=Species,y=VarP,fill=Variable))+
    geom_bar(stat="identity",position = position_stack(reverse = TRUE))+
    labs(y="Variance Proportion",fill="Predictor
(Mean Variance Proportion)",
         title="a) Presence-Absence Model")+
    scale_x_discrete(labels=Sp.labels)+
    scale_y_continuous(breaks=c(0,0.5, 1))+
    scale_fill_manual(values=customcolors,labels=var.labels)+
    theme(axis.text.x = element_text(angle=90,vjust=0),
          axis.title.x = element_blank())

  
  # Variance partitioning plot model 2 ##############
  
  m = get(modelnames[2])
  
  # grouping variables (adaptado a cada modelo)
  var.group<-attributes(m$X)$assign+1
  var.group.names<-c("Intercept","Altitude", "Lithology","Type","Age","Ext.SB","Ext.NV","Legacy")
  
  # compute model with grouping variables
  VP = computeVariancePartitioning(m,group= var.group,groupnames = var.group.names) # agrupando variables
  vals = VP$vals
  # compute explanatory powers
  preds = computePredictedValues(m)
  MF = evaluateModelFit(hM=m, predY=preds)
  
  # construyo dataframe para ggplot
  df<-melt(vals)
  names(df)<-c("Variable","Species","VarP")
  df<-droplevels(df[df$Variable!="Intercept",]) # saco el intercept
  
  #eleccion del colorpalette
  #display.brewer.all(colorblindFriendly = TRUE)
  #display.brewer.pal(n=6,name="RdBu")
  #display.brewer.pal(n=2,name="Greys")
  customcolors<-c(brewer.pal(n=6,name="RdBu"),brewer.pal(n=3,name="Greys")[1:2]) # custom colors (8) a partir de brewer colorblindfriendly
  
  # Etiquetas de las variables
  levels(df$Variable) <- c("Altitude","Lithology",
                           "Transect Type","Time since Abandonment",
                           "External TE","External NV","Legacy",
                           "Random Spatial: Plot") 
  var.vals<-ddply(df,.(Variable),summarise,val=round(mean(VarP),digits=2)) # calculo R2 promedio
  var.labels<-apply(var.vals, 1 , paste , collapse = " (") # pego el nombre de la variable con su R2
  var.labels<-paste(var.labels,")",sep="") # agrego parentesis que falta
  
  # Etiquetas de las especies
  
  fitSp<-data.frame(Species=levels(df$Species),
                    fit=round(MF$R2,digits=2)) # R2 para cada espeie
  Sp.labels<-apply(fitSp, 1 , paste , collapse = " (") # pego el nombre de la especie con su R2
  Sp.labels<-paste(Sp.labels,")",sep="") # agrego parentesis que falta
  
  # gráfico
  p2<-ggplot(data=df,aes(x=Species,y=VarP,fill=Variable))+
    geom_bar(stat="identity",position = position_stack(reverse = TRUE))+
    labs(y="Variance Proportion",fill="Predictor
(Mean Variance Proportion)",
         title="b) Abundance-COP Model")+
    scale_x_discrete(labels=Sp.labels)+
    scale_y_continuous(breaks=c(0,0.5, 1))+
    scale_fill_manual(values=customcolors,labels=var.labels)+
    theme(axis.text.x = element_text(angle=90,vjust=0))
  
  

# Variance partitioning export ####
  
p.fit<-grid.arrange(p1,p2,ncol=1)

ggsave(plot=p.fit,"./fig/fig_fit.tiff",
       width=180,height=180,units="mm",
         dpi = 600,compression="lzw")
#ggsave(plot=p.fit,"./fig/fig_fit.eps",width=180,height=180,units="mm",dpi = 600)
ggsave(plot=p.fit,"./fig/fig_fit.pdf",width=180,height=180,units="mm",
       dpi = 600, colormodel = "cmyk")

  

# Beta plots #######

# beta plot general

# PA model
  tiff("./fig/figA4.1.tif",width=180,height=140,units="mm",
     res=600,compression = "lzw") # abro el archivo

  m = get(modelnames[1])
  m$covNames <-c("(Intercept)","Altitude","Lithology_Alluvial",
                 "Lithology_Marl","Lithology_Metamorphic","Type_TL","Type_TE","Type_NV",
                 "Time","Ext.TE","Ext.NV","Legacy")
  postBeta = getPostEstimate(m, parName="Beta")
  vars = c(1:12) # selecciono las vars en m$covNames
  par(mar=c(9.1, 5.1, 2.1, 0),xaxt="n")
  plotBeta(m, post=postBeta, supportLevel = 0.95,
           param="Support",
           plotTree = FALSE,
           #covOrder = "Vector", covVector = vars,
           covNamesNumbers = c(TRUE,FALSE),
           spNamesNumbers=c(TRUE,FALSE),
           cex=c(0.9,0.9,0.7))
  title(main="Beta Plot Presence-Absence Model", line=0.5, cex.main=1.2)

  dev.off()
  
# Ab-COP model
  tiff("./fig/figA4.2.tif",width=180,height=140,units="mm",
       res=600,compression = "lzw") # abro el archivo
  
  m = get(modelnames[2])
  m$covNames <-c("(Intercept)","Altitude","Lithology_Alluvial",
                 "Lithology_Marl","Lithology_Metamorphic","Type_TL","Type_TE","Type_NV",
                 "Time","Ext.TE","Ext.NV","Legacy")
  postBeta = getPostEstimate(m, parName="Beta")
  vars = c(1:12) # selecciono las vars en m$covNames
  par(mar=c(9.1, 5.1, 2.1, 0),xaxt="n")
  plotBeta(m, post=postBeta, supportLevel = 0.95,param="Support",
           plotTree = FALSE,
           #covOrder = "Vector", covVector = vars,
           covNamesNumbers = c(TRUE,FALSE),
           spNamesNumbers=c(TRUE,FALSE),
           cex=c(0.9,0.9,0.7))
  title(main="Beta Plot Abundance-COP Model", line=0.5, cex.main=1.2)
  
  dev.off()

# para beta plot manejo local y regional por separado seleccionar las variables en cada caso
# manejo local son vars = c(7,8,9,10,11,12) # selecciono las variables en m$covNames

# Total abundance and Species richness ##############
# Bioclimates - Lithologies - Transect Type

list.bioclim<-list()
list.lito<-list()
list.type<-list()

for(i in 1:nm){
# total abundance
m = get(modelnames[i])
Gradient = constructGradient(m,focalVariable = "Altitude")
predY = predict(m, Gradient=Gradient, expected = TRUE,
                predictEtaMean = TRUE)

predS<-llply(predY,rowSums) # aplica función rowSums a cada elemento de la lista
pred<-unlist(predS,use.names = T)
pred<-data.frame(Altitude=Gradient$XDataNew$Altitude,
                 predS=as.numeric(pred))
# bioclimate
#pred<-ddply(pred,.(Bioclimate),summarise,mean=mean(predS),
#            q2.5=quantile(predS,0.025),
#            q97.5=quantile(predS,0.975))

#list.bioclim[[modelnames[i]]]<-pred

print(plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE,
                   main = modelnames[i]))

# lithology
Gradient = constructGradient(m,focalVariable = "Lithology")
predY = predict(m, Gradient=Gradient, expected = TRUE,predictEtaMean = TRUE)

predS<-llply(predY,rowSums) # aplica función rowSums a cada elemento de la lista
pred<-unlist(predS,use.names = T)
pred<-data.frame(Lithology=Gradient$XDataNew$Lithology,
                 predS=as.numeric(pred))
pred<-ddply(pred,.(Lithology),summarise,mean=mean(predS),
            q2.5=quantile(predS,0.025),
            q97.5=quantile(predS,0.975))

list.lito[[modelnames[i]]]<-pred

print(plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE,
                   main = modelnames[i]))

# transect type
Gradient = constructGradient(m,focalVariable = "Type")
predY = predict(m, Gradient=Gradient, expected = TRUE,predictEtaMean = TRUE)

predS<-llply(predY,rowSums) # aplica función rowSums a cada elemento de la lista
pred<-unlist(predS,use.names = T)
pred<-data.frame(Type=Gradient$XDataNew$Type,
                 predS=as.numeric(pred))
pred<-ddply(pred,.(Type),summarise,mean=mean(predS),
            q2.5=quantile(predS,0.025),
            q97.5=quantile(predS,0.975))

list.type[[modelnames[i]]]<-pred

print(plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE,
                   main = modelnames[i]))

}

write.csv(do.call(rbind,list.bioclim),"./panels/tab_bioclim.csv")
write.csv(do.call(rbind,list.lito),"./panels/tab_lito.csv")
write.csv(do.call(rbind,list.type),"./panels/tab_type.csv")

# prevalence - Species profiles ###################################

# Bioclimates

m = get(modelnames[1])

P.bioclimates<-data.frame(
Species=as.factor(colnames(m$Y)),
Supra=colMeans(m$Y[m$XData$Bioclimate=="Supra",]),
Meso=colMeans(m$Y[m$XData$Bioclimate=="Meso",]),
Termo=colMeans(m$Y[m$XData$Bioclimate=="Termo",])
)

P.bioclimates<-melt(P.bioclimates,value.name = "Prevalence",variable.name = "Bioclimate")

#plot grey
ggplot(data=P.bioclimates,aes(x=Species,y=Prevalence,
                              group=Bioclimate,linetype=Bioclimate))+
  geom_point()+
  geom_line()+
  scale_x_discrete(limits=names(sort(colMeans(m$Y),decreasing=T)))+
  scale_linetype_manual(values=c("solid","dashed","dotted"))+
  theme(axis.text.x = element_text(angle=90))

# plot color
p.prev1<-ggplot(data=P.bioclimates,aes(x=Species,y=Prevalence,
                              group=Bioclimate,colour=Bioclimate))+
  geom_point()+
  geom_line()+
  labs(title="a)")+
  scale_x_discrete(limits=names(sort(colMeans(m$Y),decreasing=T)))+
  theme(axis.text.x = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank())

# Lithology

m = get(modelnames[1])

P.lito<-data.frame(
  Species=as.factor(colnames(m$Y)),
  Limestone=colMeans(m$Y[m$XData$Lithology=="Limestone",]),
  Alluvial=colMeans(m$Y[m$XData$Lithology=="Alluvial",]),
  Marl=colMeans(m$Y[m$XData$Lithology=="Marl",]),
  Metamorphic=colMeans(m$Y[m$XData$Lithology=="Metamorphic",])
)

P.lito<-melt(P.lito,value.name = "Prevalence",variable.name = "Lithology")

#plot grey
p.prev2<-ggplot(data=P.lito,aes(x=Species,y=Prevalence,
                              group=Lithology,linetype=Lithology))+
  geom_point()+
  geom_line()+
  labs(title="b)")+
  scale_x_discrete(limits=names(sort(colMeans(m$Y),decreasing=T)))+
  scale_linetype_manual(values=c("solid","dashed","dotdash","dotted"))+
  theme(axis.text.x = element_text(angle=90),legend.position="top")

# plot color
ggplot(data=P.lito,aes(x=Species,y=Prevalence,
                              group=Lithology,colour=Lithology))+
  geom_point()+
  geom_line()+
  scale_x_discrete(limits=names(sort(colMeans(m$Y),decreasing=T)))+
  theme(axis.text.x = element_text(angle=90),legend.position="top")

#export
p.prev<-grid.arrange(p.prev1,p.prev2,heights = c(0.45, 0.55))

ggsave(plot=p.prev,"./fig/fig_prevalence.tiff",width=180,height=140,units="mm",
       dpi = 600,compression="lzw")
ggsave(plot=p.prev,"./fig/fig_prevalence.eps",width=180,height=140,units="mm",dpi = 600)
ggsave(plot=p.prev,"./fig/fig_prevalence.pdf",width=180,height=140,units="mm",
       dpi = 600, colormodel = "cmyk")



# Age predictions ###############################################

for(j in 1:nm){
m = get(modelnames[j])

Gradient = constructGradient(m,focalVariable = "Age")
predY = predict(m, Gradient=Gradient, expected = TRUE,predictEtaMean = TRUE)

predS<-llply(predY,rowSums) # aplica función rowSums a cada elemento de la lista
pred<-unlist(predS,use.names = T)
pred<-data.frame(Age=Gradient$XDataNew$Age,
                 predS=as.numeric(pred))
pred<-ddply(pred,.(Age),summarise,mean=mean(predS),
            q2.5=quantile(predS,0.025),
            q97.5=quantile(predS,0.975))
pred.general<-pred

if(j==1){
  pS<-ggplot()+
  geom_line(data=pred.general, aes(x=Age, y=mean),size=1.4)+
  geom_ribbon(data=pred.general, aes(x=Age, ymin=q2.5,ymax=q97.5),
              alpha=0.1)+
  labs(title="a)",y="Species Richness")+
  theme(axis.title.x = element_blank())
  #geom_line(data=pred.bioclimate, aes(x=Age, y=mean,linetype=Bioclimate),colour="darkgrey")+
  #scale_linetype_manual(values=c("dotdash", "dashed", "dotted"))
}else{
pAb<-ggplot()+
  geom_line(data=pred.general, aes(x=Age, y=mean),size=1.4)+
  geom_ribbon(data=pred.general, aes(x=Age, ymin=q2.5,ymax=q97.5),
              alpha=0.1)+
  labs(title="b)",y="Total Abundance (ind. / transect)",x="Time since Abandonment (Years)")
}
}

p_age<-grid.arrange(pS,pAb)

# export
ggsave(plot=p_age,"./fig/fig_age.tiff",width=180,height=140,units="mm",
       dpi = 600,compression="lzw")

ggsave(plot=p_age,"./fig/fig_age.eps",width=180,height=140,units="mm",dpi = 600)

ggsave(plot=p_age,"./fig/fig_age.pdf",width=180,height=140,units="mm",dpi = 600, 
       colormodel = "cmyk")


# Predictions Tree Lines - Terrace - Natural Patch #########################################

# VERSION DE PRUEBA TREELINE

m = get(modelnames[1])
postBeta = getPostEstimate(m, parName="Beta")
Gradient <- constructGradient(m,focalVariable = "Type",
                              non.focalVariables=list("Altitude"=list(2),
                                                      "Lithology"=list(2)))

Gradient <- constructGradient(m,focalVariable = "Type",
                              non.focalVariables=list("Altitude"=list(3,1200),
                                                      "Lithology"=list(1)))
summary(m$XData$Altitude)
#"Ext.SB"=list(2)

#Gradient$XDataNew<-Gradient$XDataNew[1:2,]
#Gradient$studyDesignNew<-Gradient$studyDesignNew[1:2,]

predY = predict(m, Gradient=Gradient, expected = TRUE,
                predictEtaMean = TRUE)  

pred<-melt(predY)[1:3]
names(pred)<-c("Type","Sp","pred")
pred<-pred[pred$Type%in%1:2,]
pred$Type<-factor(pred$Type,labels = c("Open Space","Tree Line"))


pred<-ddply(pred,.(Type,Sp),summarise,mean=mean(pred),
            q2.5=quantile(pred,0.025),
            q97.5=quantile(pred,0.975))

selPar<-which(m$covNames=="TypeTL") # seleccionar variable para asteriscos
sig<-postBeta$support[selPar,]>0.95 # overall beta significativo (para los tres bioclims juntos)
sig<-factor(sig,labels=c("","*")) # armo las labels para el grafico

ggplot(data=pred,aes(x=Sp,colour=Type))+
  geom_pointrange(aes(y=mean,min=q2.5,max=q97.5),position = position_dodge(width=0.2))+
  annotate(geom="text",x=names(sig),y = 1.0,label = sig,size=6)+
  labs(y="Probability of Occurrence",
       x="Species")+
  scale_colour_manual(values=c("grey","black"))+
  theme(axis.text.x = element_text(angle=90),
        legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0,5.5,5.5,5.5), "points"))

# VERSION CON TODAS

m = get(modelnames[1])
postBeta = getPostEstimate(m, parName="Beta")


    #Gradient <- constructGradient(m,focalVariable = "Type",
    #                            non.focalVariables=list("Bioclimate"=list(3,b),
    #                                                    "Lithology"=list(1)))
    
    Gradient <- constructGradient(m,focalVariable = "Type",
                                  non.focalVariables=list("Altitude"=list(2),
                                                          "Lithology"=list(2)))
                                                          #"Ext.SB"=list(2),
                                                          #"Ext.NV"=list(2)))
    
    Gradient2 <- constructGradient(m,focalVariable = "Type",
                                   non.focalVariables=list("Altitude"=list(3,1200),
                                                           "Lithology"=list(1)))
                                                           #"Ext.SB"=list(2),
                                                           #"Ext.NV"=list(2)))
    # predicciones sobre Limestone para la mayoría de las especies
    #Gradient$XDataNew<-Gradient$XDataNew[1:2,]
    #Gradient$studyDesignNew<-Gradient$studyDesignNew[1:2,]
    predY = predict(m, Gradient=Gradient, expected = TRUE,
                    predictEtaMean = TRUE)  
    pred<-melt(predY)[1:3]
    names(pred)<-c("Type","Sp","pred")
    pred<-pred[pred$Type%in%1:2,]
    pred$Type<-factor(pred$Type,labels = c("Open Space","Tree Line"))
    
    # predicciones sobre Metamorphic para AspAcu y AspAlb
    #Gradient2$XDataNew<-Gradient2$XDataNew[1:2,]
    #Gradient2$studyDesignNew<-Gradient2$studyDesignNew[1:2,]
    predY2 = predict(m, Gradient=Gradient2, expected = TRUE,
                     predictEtaMean = TRUE)  
    pred2<-melt(predY2)[1:3]
    names(pred2)<-c("Type","Sp","pred")
    pred2<-pred2[pred2$Type%in%1:2,]
    pred2$Type<-factor(pred2$Type,labels = c("Open Space","Tree Line"))
    
    # junto las planillas
    
    #pred2<-pred2[pred2$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
    #pred<-pred[!pred$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
    pred2<-pred2[pred2$Sp%in%c("DapGni","JunOxy","JunPho","RosSp"),]
    pred<-pred[!pred$Sp%in%c("DapGni","JunOxy","JunPho","RosSp"),]
    
    pred<-rbind(pred,pred2)
    
    pred<-ddply(pred,.(Type,Sp),summarise,mean=mean(pred),
                q2.5=quantile(pred,0.025),
                q97.5=quantile(pred,0.975))
    
    selPar<-which(m$covNames=="TypeTL") # seleccionar variable para asteriscos
    sig<-postBeta$support[selPar,]>0.95 # overall beta significativo (para los tres bioclims juntos)
    sig<-factor(sig,labels=c("","*")) # armo las labels para el grafico

  p_TL <- ggplot(data=pred,aes(x=Sp,colour=Type))+
      geom_pointrange(aes(y=mean,min=q2.5,max=q97.5),position = position_dodge(width=0.2))+
      annotate(geom="text",x=names(sig),y = 0.99,label = sig,size=6)+
      labs(subtitle="a)",
           y="",
           x="")+
    scale_y_continuous(limits=c(0,1.02))+
      scale_colour_manual(values=c("grey","black"))+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            plot.margin = unit(c(0,5.5,5.5,5.5), "points"))


# External Source - Terrace Edge ###########################
# VERSION POR SEPARADO

  Gradient <- constructGradient(m,focalVariable = "Ext.SB",
                                non.focalVariables=list("Altitude"=list(3,1100),
                                                        "Lithology"=list(2)))
                                                        #"Ext.NV"=list(2)))
                                                        #"Type"=list(3,"OS")))
  
  Gradient2 <- constructGradient(m,focalVariable = "Ext.SB",
                                 non.focalVariables=list("Altitude"=list(2),
                                                         "Lithology"=list(3,"Metamorphic")))
                                                         #"Ext.NV"=list(2)))
                                                         #"Type"=list(3,"OS")))
  # predicciones sobre Limestone para la mayoría de las especies
  #Gradient$XDataNew<-Gradient$XDataNew[1:2,]
  #Gradient$studyDesignNew<-Gradient$studyDesignNew[1:2,]
  predY = predict(m, Gradient=Gradient, expected = TRUE,
                  predictEtaMean = TRUE)  
  pred<-melt(predY)[1:3]
  names(pred)<-c("Type","Sp","pred")
  pred<-pred[pred$Type%in%1:2,]
  pred$Type<-factor(pred$Type,labels = c("Without Terrace Edge","With Terrace Edge"))
  
  # predicciones sobre Metamorphic para AspAcu y AspAlb
  #Gradient2$XDataNew<-Gradient2$XDataNew[1:2,]
  #Gradient2$studyDesignNew<-Gradient2$studyDesignNew[1:2,]
  predY2 = predict(m, Gradient=Gradient2, expected = TRUE,
                   predictEtaMean = TRUE)  
  pred2<-melt(predY2)[1:3]
  names(pred2)<-c("Type","Sp","pred")
  pred2<-pred2[pred2$Type%in%1:2,]
  pred2$Type<-factor(pred2$Type,labels = c("Without Terrace Edge","With Terrace Edge"))
  
  # junto las planillas
  
  #pred2<-pred2[pred2$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
  #pred<-pred[!pred$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
  pred2<-pred2[pred2$Sp%in%c("AspAlb","AspHor"),]
  pred<-pred[!pred$Sp%in%c("AspAlb","AspHor"),]
  
  pred<-rbind(pred,pred2)
  
  pred<-ddply(pred,.(Type,Sp),summarise,mean=mean(pred),
              q2.5=quantile(pred,0.025),
              q97.5=quantile(pred,0.975))
  
  selPar<-which(m$covNames=="Ext.SB1") # seleccionar variable para asteriscos
  sig<-postBeta$support[selPar,]>0.95 # overall beta significativo (para los tres bioclims juntos)
  sig<-factor(sig,labels=c("","*")) # armo las labels para el grafico

  p_TE <- ggplot(data=pred,aes(x=Sp,colour=Type))+
    geom_pointrange(aes(y=mean,min=q2.5,max=q97.5),position = position_dodge(width=0.2))+
    annotate(geom="text",x=names(sig),y = 0.99,label = sig,size=6)+
    labs(subtitle="b)",
         y="Probability of Occurrence",
         x="")+
  scale_y_continuous(limits=c(0,1.02))+
  scale_colour_manual(values=c("grey","black"))+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          plot.margin = unit(c(0,5.5,5.5,5.5), "points"))

# External Source - Natural Vegetation ###########################
# VERSION POR SEPARADO

  Gradient <- constructGradient(m,focalVariable = "Ext.NV",
                                non.focalVariables=list("Altitude"=list(2),
                                                        "Lithology"=list(2)))
                                                        #"Type"=list(3,"OS")))
  
  #Gradient2 <- constructGradient(m,focalVariable = "Ext.NV",
  #                               non.focalVariables=list("Bioclimate"=list(3,b),
  #                                                       "Lithology"=list(3,"Metamorphic")))
                                                         #"Type"=list(3,"OS")))

  # predicciones sobre Limestone para la mayoría de las especies
  #Gradient$XDataNew<-Gradient$XDataNew[1:2,]
  #Gradient$studyDesignNew<-Gradient$studyDesignNew[1:2,]
  predY = predict(m, Gradient=Gradient, expected = TRUE,
                  predictEtaMean = TRUE)  
  pred<-melt(predY)[1:3]
  names(pred)<-c("Type","Sp","pred")
  pred<-pred[pred$Type%in%1:2,]
  pred$Type<-factor(pred$Type,labels = c("Without Natural Vegetation Patch","With Natural Vegetation Patch"))
  
  # predicciones sobre Metamorphic para AspAcu y AspAlb
  #Gradient2$XDataNew<-Gradient2$XDataNew[1:2,]
  #Gradient2$studyDesignNew<-Gradient2$studyDesignNew[1:2,]
  #predY2 = predict(m, Gradient=Gradient2, expected = TRUE,
  #                 predictEtaMean = TRUE)  
  #pred2<-melt(predY2)[1:3]
  #names(pred2)<-c("Type","Sp","pred")
  #pred2$Type<-factor(pred2$Type,labels = c("Without Natural Vegetation Patch","With Natural Vegetation Patch"))
  
  # junto las planillas
  
  #pred2<-pred2[pred2$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
  #pred<-pred[!pred$Sp%in%c("AspAcu","AspAlb"),] # antes usaba metamorfic para ambas
  #pred2<-pred2[pred2$Sp%in%c("AspAcu"),]
  #pred<-pred[!pred$Sp%in%c("AspAcu"),]
  #pred<-rbind(pred,pred2)
  
  pred<-ddply(pred,.(Type,Sp),summarise,mean=mean(pred),
              q2.5=quantile(pred,0.025),
              q97.5=quantile(pred,0.975))
  
  selPar<-which(m$covNames=="Ext.NV1") # seleccionar variable para asteriscos
  sig<-postBeta$support[selPar,]>0.95 # overall beta significativo (para los tres bioclims juntos)
  sig<-factor(sig,labels=c("","*")) # armo las labels para el grafico

  p_NV <- ggplot(data=pred,aes(x=Sp,colour=Type))+
    geom_pointrange(aes(y=mean,min=q2.5,max=q97.5),position = position_dodge(width=0.2))+
    annotate(geom="text",x=names(sig),y = 0.99,label = sig,size=6)+
    labs(subtitle="c)",
         y="",
         x="Species")+
    scale_colour_manual(values=c("grey","black"))+
    scale_y_continuous(limits=c(0,1.02))+
    theme(axis.text.x = element_text(angle=90),
          legend.position = "top",
          legend.title = element_blank(),
          plot.margin = unit(c(0,5.5,5.5,5.5), "points"))

p_pred <- grid.arrange(p_TL,p_TE,p_NV,
             heights=c(0.3,0.3,0.36))
             #widths=c(0.35,0.35,0.3))
  
#export
ggsave(plot=p_pred,"./fig/fig_TLTENV.tiff",width=180,height=210,units="mm",dpi = 600,compression="lzw")
ggsave(plot=p_pred,"./fig/fig_TLTENV.eps",width=180,height=210,units="mm",dpi = 600)
ggsave(plot=p_pred,"./fig/fig_TLTENV.pdf",width=180,height=210,units="mm",colormodel = "cmyk",dpi = 600)



### distribución altitudinal de cada especie ########################

m = get(modelnames[1])
Gradient = constructGradient(m,focalVariable = "Altitude")
predY = predict(m, Gradient=Gradient, expected = TRUE)  

    pred<-melt(predY)[1:3]
    names(pred)<-c("Altitude","Sp","pred")
    
    pred<-ddply(pred,.(Altitude,Sp),summarise,mean=mean(pred),
                q2.5=quantile(pred,0.025),
                q97.5=quantile(pred,0.975))
    
    p_sp_alt<-ggplot()+
        geom_line(data=pred, aes(x=Altitude, y=mean),size=1.4)+
        geom_ribbon(data=pred, aes(x=Altitude, ymin=q2.5,ymax=q97.5),
                    alpha=0.1)+
        labs(y="Probability of Ocurrence")+
      facet_wrap(~Sp,ncol=3)
    
    #export
    ggsave(plot=p_sp_alt,"./fig/fig_sp_altitude.tiff",width=180,height=140,units="mm",dpi = 600,compression="lzw")
    ggsave(plot=p_sp_alt,"./fig/fig_sp_altitude.eps",width=180,height=140,units="mm",dpi = 600)
    ggsave(plot=p_sp_alt,"./fig/fig_sp_altitude.pdf",width=180,height=140,units="mm",colormodel = "cmyk",dpi = 600)
    
    
### Biologial traits ###############################################

m = get(modelnames[1])

# traits vs altitude
Gradient = constructGradient(m,focalVariable = "Altitude")
predY = predict(m, Gradient=Gradient, expected = TRUE)  

colnames(m$Tr)

ylab <- c("Growth Form Liana (Proportion)",
          "Growth Form Shrub (Proportion)",
          "Root Depth (m)",
          "Mixed Dispersal (Proportion)",
          "Fruit length (mm)",
          "Fruit Water Content (%)")

tiff("./fig/fig_traits_altitude.tif",width=180,height=140,units="mm",
     res=600,compression = "lzw") # abro el archivo


par(mfrow=c(3,2))
j=1

for (i in c(2,3,4,7,8,10)){
  plotGradient(m, Gradient, pred=predY, measure="T",
                    index=i,
               ylabel=ylab[j],
               xlabel="Altitude (m a.s.l.)",
             #showData = TRUE,
             #yshow = 0,
             #main = ": community weighted mean trait (total effect)",
             showPosteriorSupport = FALSE)
  j <- j+1
}
dev.off()
# guardar plot manualmente


# traits vs age
Gradient = constructGradient(m,focalVariable = "Age")
predY = predict(m, Gradient=Gradient, expected = TRUE)  

colnames(m$Tr)

ylab <- c("Root Depth (m)",
          "Mixed Dispersal (Proportion)",
          "Fruit length (mm)",
          "Seed Dry Mass (g)",
          "Fruit Water Content (%)")

tiff("./fig/fig_traits_age.tif",width=180,height=140,units="mm",
     res=600,compression = "lzw") # abro el archivo

par(mfrow=c(3,2))
j=1
for (i in c(4,7,8,9,10)){
  plotGradient(m, Gradient, pred=predY, measure="T",
               index=i,
               ylabel=ylab[j],
               xlabel="Time Since Abandonment (years)",
               #showData = TRUE,
               #yshow = 0,
               #main = ": community weighted mean trait (total effect)",
               showPosteriorSupport = FALSE)
  j <- j+1
}
# guardar plot manualmente
dev.off()

# Gamma Plot ####

# PA model
tiff("./fig/figA_gammaplotPA.tif",width=180,height=140,units="mm",
     res=600,compression = "lzw") # abro el archivo

m = get(modelnames[1])

m$trNames<-c("(Intercept)","Growth Form Liana",
  "Growth Form Shrub",
  "Root Depth","Specific Leaf Area",
  "Sexual System Monoecious",
  "Mixed Dispersal",
  "Fruit length",
  "Seed Dry Mass",
  "Fruit Water Content")

m$covNames <-c("(Intercept)","Altitude","Lithology_Alluvial",
               "Lithology_Marl","Lithology_Metamorphic","Type_TL","Type_TE","Type_NV",
               "Time","Ext.TE","Ext.NV","Legacy")

postGamma = getPostEstimate(m, parName="Gamma")
plotGamma(m, post=postGamma, supportLevel = 0.95, param="Sign",
          covNamesNumbers = c(TRUE,FALSE),
          #trNamesNumbers=c(m$nt<21,FALSE),
          trNamesNumbers=c(TRUE,FALSE),
          cex=c(0.6,0.6,0.8))
title(main="Gamma Plot Presence-Absence Model", line=2.5,cex.main=0.8)


dev.off()
