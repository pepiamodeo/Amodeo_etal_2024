# Data and R scripts for the article Amodeo *et al* 2024 entitled *Bridging the gap between ecological succession of fleshy-fruited shrubs and restoration frameworks in semiarid oldfields*.

Authors: 

* Amodeo, M.
* Martínez-López, V.
* Zapata-Pérez, V.
* Robledano-Aymerich, F. 

This repository is preserved by Zenodo. [![DOI](https://zenodo.org/badge/768723113.svg)](https://zenodo.org/doi/10.5281/zenodo.13127173)

### Project description

This R project consists in the datasets and scripts used for the analyses and plots exposed in the article. 

Data were analyzed using Hierarchical Modelling of Species Communities (HMSC Ovaskainen & Abrego, 2020), a multivariate hierarchical generalized linear mixed model fitted with Bayesian inference. Data comprise the abundances of fleshy-fruited shrub species surveyed in 299 transects distributed over 116 oldfields, each field hosting one to four transects. We included information about functional and life-history traits in a joint species modeling approach to analyze how species traits affect the natural recolonization of shrubs in oldfields. The count of each of the species was used as the response variable and, because of the zero-inflation nature of the data, we applied a Hurdle model approach: one model (probit regression approach) for the occurrence data (presence-absence) and another (normal linear regression approach) for the abundance conditional on presence (henceforth abundance-COP model), declaring zeros as missing data, log-transforming and scaling to zero mean and unit variance within each species. Climatic numerical variables obtained for each field from the EuMedClim database were included as predictor variables (Fréjaville & Benito Garzón, 2018). We fitted the HMSC model using the R-package Hmsc (Tikhonov et al., 2020). Custom plots were constructed using the package ggplot2 (Wickham, 2009) in R (R Core Team, 2020). 

For a full description of the study, data, experimental design, analyses and conclusions please see [DOI].


# References

* Fréjaville, T., Benito Garzón, M., 2018. The EuMedClim Database: Yearly Climate Data (1901–2014) of 1 km Resolution Grids for Europe and the Mediterranean Basin. Front. Ecol. Evol. 6, 31. https://doi.org/10.3389/fevo.2018.00031

* Ovaskainen, O., Abrego, N., 2020. Joint Species Distribution Modelling: With Applications in R, Ecology, Biodiversity and Conservation. Cambridge University Press, Cambridge. https://doi.org/10.1017/9781108591720

* R Core Team, 2020. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.

* Tikhonov, G., Opedal, Ø.H., Abrego, N., Lehikoinen, A., Jonge, M.M.J., Oksanen, J., Ovaskainen, O., 2020. Joint species distribution modelling with the r ‐package H msc. Methods in Ecology and Evolution 11, 442–447. https://doi.org/10.1111/2041-210X.13345

* Wickham, H., 2009. ggplot2: Elegant Graphics for Data Analysis. Springer New York, New York, NY. https://doi.org/10.1007/978-0-387-98141-3
