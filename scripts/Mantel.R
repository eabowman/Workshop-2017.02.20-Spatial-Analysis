## Script created by Liz Bowman, 2017-02-15 for use in a workshop on how to examine
## data for spatial paterns and to test for correlations between communities
## and environmental variables

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
# Mantel test
#=========================================================================================

#-----------------------------------------------------------------------------------------
# Check for spatial autocorrelation in mite community data
#-----------------------------------------------------------------------------------------

#<< Create distance matrices >> ----------------------------------------------------------
#--geographical distance matrix using Euclidean distance default
mite.dist.spatial <- dist(mite.xy)

#--community distance matrix using Morisitia's index
mite.dist.comm <- vegdist(mite, method = 'morisita')

#<< Run mantel test >> -------------------------------------------------------------------
#--Null hypothesis: community similarity is indepent of geographical proximity
#--Alternative hypothesis: communities closer geographically are more similar
mite.spatial <- mantel(mite.dist.spatial, mite.dist.comm)
mite.spatial

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot of community as a function of geographic distance classes
mite.spatial.corr <- mantel.correlog(mite.dist.comm, mite.dist.comm)
plot(mite.spatial.corr)

#--determine where the break points
mite.spatial.corr$break.pts

#-----------------------------------------------------------------------------------------
# Check for correlation between multiple environmental variables and community data
#-----------------------------------------------------------------------------------------
#<< Create distance matrices >> ----------------------------------------------------------
#--community distance matrix
mite.dist.comm <- vegdist(mite, method = 'morisita')

#--Environmental distance matrix

#--Substrate: a categorical variable
unique(mite.env$Substrate)
#--we need to change this over to numbers that represent each level
within(mite.env, Substrate <- factor(Substrate, labels = c(1:7)))
mite.env$Substrate <- as.numeric(mite.env$Substrate)

mite.dist.env <- dist(mite.env[1:3])

#<< Run mantel test of environment and community data >> ----------------------------------
mite.env.comm <- mantel(mite.dist.env, mite.dist.comm)
mite.env.comm

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot community dissimilarityy as a function of environment
mite.env.comm.corr <- mantel.correlog(mite.dist.env, mite.dist.comm)
plot(mite.env.comm.corr)

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
#<< Create distance matrices >> ----------------------------------------------------------
#--ozone distance matrix
mite.dist.comm <- vegdist(mite, method = 'morisita')

#--Water content of the substrate (g/L): a continuous variable
mite.dist.watrcont <- dist(mite.env$WatrCont)

#<< Run mantel test >> -------------------------------------------------------------------
mite.watr.comm <- mantel(mite.dist.watrcont, mite.dist.comm)
mite.watr.comm

#--compare to community data and mult. environmental variables
mite.env.comm

#<< Mantel correlogram >> ----------------------------------------------------------------
#--correlograms plot of community as a function of geographic distance classes
mite.watrcont.comm.corr <- mantel.correlog(mite.dist.watrcont, mite.dist.comm)
plot(mite.watrcont.comm.corr)

#--determine where the break points
mite.watrcont.comm.corr$break.pts
