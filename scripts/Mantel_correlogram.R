## Script created by Liz Bowman, 2017-02-15 for use in a workshop on how to examine
## data for spatial paterns and to test for correlations between communities
## and environmental variables
## eabowman@email.arizona.edu

#=========================================================================================
# Load packages and data
#=========================================================================================

#<< packages for analyses >> -------------------------------------------------------------
#--Uncomment if you have not previously installed these programs
#install.packages('vegan')

#--once above packages are installed run the following code to load them
library(vegan)

#-----------------------------------------------------------------------------------------
# input data for analyses
#-----------------------------------------------------------------------------------------

#<< Mite data >> -------------------------------------------------------------------------
#--community data
data(mite)
#--environmental data
data(mite.env)
#--location data
data(mite.xy)

#--if you are interested in what all of the data is then check out the help page
#?mite

#=========================================================================================
# Mantel correlogram
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Check for spatial autocorrelation in mite community data
#-----------------------------------------------------------------------------------------

#<< Create distance matrices >> ----------------------------------------------------------
#--geographical distance matrix using Euclidean distance default
mite.dist.spatial <- dist(mite.xy)

#--detrend the species data by regresion on the site coordinates
mite.hel <- decostand(mite, 'hellinger')
#--test for linear trends; if signficiant, then we will need to detrend the data
anova(rda(mite.h,mite.xy))

#--to detrend, we will regress all variables against the x-y coordinates and retain
#--the residuals
#--computation of linearly detrended mite data
mite.hel.resid <- resid(lm(as.matrix(mite.hel) ~ ., data = mite.xy))

#--recalculate distance matrix using detrended data
mite.hel.D <- dist(mite.hel.resid)

#<< Run mantel test >> -------------------------------------------------------------------
#--Null hypothesis: community similarity is indepent of geographical proximity
#--Alternative hypothesis: communities closer geographically are more similar
mite.spatial <- mantel(mite.dist.spatial, mite.hel.D)
mite.spatial

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot of community as a function of geographic distance classes
# if you want to specify break points use the argument break.pts and supply a vector with
# the break pts you want the function to use
# cutoff = T which is the default limts the correlogram to the distance classes that
# include all pts
mite.spatial.corr <- mantel.correlog(mite.dist.spatial, mite.hel.D, cutoff = F)
plot(mite.spatial.corr)

#--determine where the break points
mite.spatial.corr$break.pts

#--print the p-values for each distance class
mite.spatial.corr$mantel.res

#-----------------------------------------------------------------------------------------
# Check for correlation between multiple environmental variables and community data
#-----------------------------------------------------------------------------------------
#<< Create distance matrices >> ----------------------------------------------------------

### Environmental distance matrix
#--Substrate: a categorical variable
unique(mite.env$Substrate)
#--we need to change this over to numbers that represent each level
within(mite.env, Substrate <- factor(Substrate, labels = c(1:7)))
mite.env$Substrate <- as.numeric(mite.env$Substrate)
#--make distance matrix
mite.dist.env <- dist(mite.env[1:3])

#<< Run mantel test of environment and community data >> ---------------------------------
mite.env.comm <- mantel(mite.dist.env, mite.hel.D)
mite.env.comm

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot community dissimilarityy as a function of environment
mite.env.comm.corr <- mantel.correlog(mite.dist.env, mite.hel.D, cutoff = F)
plot(mite.env.comm.corr)

#--determine where the break points
mite.env.comm.corr$break.pts

#--print the p-values for each distance class
mite.env.comm.corr$mantel.res

#-----------------------------------------------------------------------------------------
# Check for correlation between single environmental variables and geographic distance
#-----------------------------------------------------------------------------------------

#<< Create distance matrices >> ----------------------------------------------------------
#--geographical distance matrix
mite.dist.spatial <- dist(mite.xy)

#--Water content of the substrate (g/L): a continuous variable
mite.dist.watrcont <- dist(mite.env$WatrCont)

#--Substrate: a categorical variable
unique(mite.env$Substrate)
#--we need to change this over to numbers that represent each level
within(mite.env, Substrate <- factor(Substrate, labels = c(1:7)))
mite.env$Substrate <- as.numeric(mite.env$Substrate)
#---now we can make a distance matrix of it
mite.dist.subs <- dist(mite.env$Substrate)

#<< Run mantel test of water content of substrate as a function of geographical distance>>
mite.watrcont.geo <- mantel(mite.dist.spatial, mite.dist.watrcont)
mite.watrcont.geo

#<< Run mantel test of water content of substrate as a function of geographical distance>>
mite.subs.geo <- mantel(mite.dist.spatial, mite.dist.subs)
mite.subs.geo

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot water content of the substrate as a function of distance
mite.watrcont.corr <- mantel.correlog(mite.dist.watrcont, mite.dist.spatial,
                                      cutoff = F) #
plot(mite.watrcont.corr)

#--determine where the break points
mite.watrcont.corr$break.pts

#--which points are significant
mite.watrcont.corr$mantel.res

#--correlograms plot water content of the substrate as a function of distance
mite.subs.corr <- mantel.correlog(mite.dist.subs, mite.dist.spatial)
plot(mite.subs.corr)

#-----------------------------------------------------------------------------------------
# Check for correlation between water content of substrate and community
#-----------------------------------------------------------------------------------------
#<< Run mantel test >> -------------------------------------------------------------------
mite.watr.comm <- mantel(mite.dist.watrcont, mite.hel.D)
mite.watr.comm

#--compare to community data and mult. environmental variables
mite.env.comm

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot of community as a function of geographic distance classes
mite.watrcont.comm.corr <- mantel.correlog(mite.dist.watrcont, mite.hel.D)
plot(mite.watrcont.comm.corr)

#--compare to community data and mult. environmental variables
# Based on the similarity in r-values and correlograms, it seems that water content of
# the substrate explains more of the variation in the community data than substrate
# density or substrate
plot(mite.env.comm.corr) 

#--determine where the break points
mite.watrcont.comm.corr$break.pts

#--which points are significant
mite.watrcont.comm.corr$mantel.res
