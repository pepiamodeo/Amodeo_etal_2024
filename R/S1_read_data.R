# MODIFY THIS SCRIPT SO THAT IT RUNS WITH YOUR OWN DATA

# You need to provide an SXY file.
# The files TP and P are optional, so indicate with TRUE/FALSE if they are included or not
is.TP = TRUE
is.P = FALSE

# READING IN SXY: study design (S) and/or covariates (X) and species data (Y)
SXY = read.csv("./data/SXY.csv", stringsAsFactors=TRUE) # or use read.csv2 when values are semicolon-separated
Amb = read.csv("./data/Amb.csv", stringsAsFactors=TRUE) # or use read.csv2 when values are semicolon-separated

SXY<-merge(SXY,Amb,by="PlotID")

#### INICIO SEPT 2022
# genero variable indicador binario de legacy: cultivo de arbol o cultivo no arbol

df<-reshape2::dcast(SXY, PlotID ~ Type,
                      fun.aggregate = length,
                      value.var="Age")[,c(1,5)] # defino presencia de TL en la parcela
names(df)[2]<-"Legacy" # la renombro legacy TL

SXY<-merge(SXY,df,by="PlotID") # incorporo ese indicador a X

# filtro solo transectos dentro del campo OS y TL, saco SB y NV
SXY_filt<-dplyr::filter(SXY,Type%in%c("OS","TL"))
SXY_filt<-droplevels(SXY_filt)

levels(SXY_filt$PlotID)
levels(SXY$PlotID)

#### FIN SEPT 2022

# Modify the next three lines to split your SXY file to components that relate to
# S: study design, including units of study and their possible coordinates (named as Route_x and Route_y to indicate that they relate to the units of Route)
# X: covariates to be used as predictors
# Y: species data
# If you don't have variables that define the study design, indicate this by S=NULL
# If you don't have covariate data, indicate this by X=NULL

names(SXY)
S=SXY[,c(1,2,35,36)] # podría agergar bioclim y litology para tener opcion de incorporarlos random
X=SXY[,c(3,4,5,6,7,8,37,38,39,40,41)]
Y=SXY[,9:34]

S_filt=SXY_filt[,c(1,2,35,36)] # podría agergar bioclim y litology para tener opcion de incorporarlos random
X_filt=SXY_filt[,c(3,4,5,6,7,8,37,38,39,40,41)]
Y_filt=SXY_filt[,9:34]

# What is not always easy is to decide what goes to S and what to X.
# As a general rule, include in S those variables that you think should be modelled as random effect,
# and in X those that you think should be modelled as fixed effects.
# Don't worry if you are not sure which one is the "right choice", we will discuss this with you.

# Check that the data looks as it should!
#View(S)
#View(X)
#View(Y)
# check that community data are numeric and have finite numbers. If the script
# writes "Y looks OK", you are ok.
if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
    print("Y looks OK")
} else {
	print("Y should be numeric and have finite values")	}
# Check that the stydy design data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(S))) {
  print("S has NA values - not allowed for")
} else {
  print("S looks ok")	}
# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(X))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}


# READING IN TP: traits (T) and/or phylogenetic information in table format (P)
if(is.TP){
  # Read in the species names as rownames, not as a column of the matrix
  TP = read.csv("./data/Tr.csv", stringsAsFactors=TRUE,row.names = 2)
  # The script below checks if the species names in TP are identical and in the same order as in Y
  # If the script prints "species names in TP and SXY match", you are ok.
  # If it says that they do not match, you need to modify the files so that they match
  if(all(rownames(TP)==colnames(Y))) {
    print("species names in TP and SXY match")
  } else{
    print("species names in TP and SXY do not match")
  }
  # Modify the next two lines to split your TP file to components that relate to
  # Tr: species traits (note that T is a reserved word in R and that's why we use Tr)
  # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species
  # If you don't have trait data, indicate this by Tr=NULL.
  # If TP does not have phylogenetic data (because you don't have such data at all, or because
  # it is given in tree-format, like is the case in this example), indicate this with P=NULL
  names(TP)
  Tr = TP[,c(8:20)]
  P = TP[,c(1,4:6)]
  # Check that the data looks as it should!
  #View(Tr)
  #View(P)
  # Check that the Tr data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(Tr))) {
    print("Tr has NA values - not allowed for")
  } else {
    print("Tr looks ok")	}
  # Check that the phylogenetic/taxonomic data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(P))) {
    print("P has NA values - not allowed for")
  } else {
    print("P looks ok")	}
}

# Chequeo NA en Tr para las 12 sp seleccionadas
is.na(Tr)[spSEL,] # traer spSEL de script S2

apply(X=Tr,MARGIN=1,FUN=anyNA)[spSEL]

# READING IN P: phylogenetic information in tree format (P)
# we use ape package for trees, and P.tre must be in a format that ape understands
if(is.P){
  # Read in the phylogenetic tree using read.tree from ape
  library(ape)
  P = read.tree("./analyses/HMSC/P.tre")
  # When you look at P (e.g. write P and press enter),
  # you should see that it is a phylogenetic tree which
  # is rooted and includes branch lengths and tip labels
  # The script below checks if the species names in P are identical (but not necessarily in the same order) as in Y
  # If the script prints "species names in P and SXY match", you are ok.
  # If it says that they do not match, you need to modify the files so that they match
  if(all(sort(P$tip.label) == sort(colnames(Y)))){
    print("species names in P and SXY match")
  } else{
    print("species names in P and SXY do not match")
  }
  # Check that the data looks as it should!
  plot(P, cex=0.5)
}
