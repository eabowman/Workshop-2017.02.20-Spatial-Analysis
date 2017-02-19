## Script created by Liz Bowman, 2017-02-15 for use in a workshop on how to examine
## data for spatial paterns and to test for correlations between communities
## and environmental variables

#=========================================================================================
# Load packages and data
#=========================================================================================

#<< packages for analyses >> -------------------------------------------------------------
#--Uncomment if you have not previously installed these programs
#install.packages('vegan')
#install.packages('spdep)

#--once above packages are installed run the following code to load them
library(vegan)
library(spdep)

#--load the function plot.links from script folder
source('./scripts/plot.links.R')

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

#=========================================================================================
# Spatial correlogram based on Moran's I
#=========================================================================================

#--normalize community data so that it can be used with euclidean distance
mite.h <- decostand(mite,'hellinger')

#--Center the spatial data
mite.xy.c <- scale(mite.xy, center=T, scale=F)

#-----------------------------------------------------------------------------------------
# Search for neighbors of all points wihtin a radius of 0.7 m and multiples (i.e., 0 to 
# 0.7m, 0.7 to 1.4m, etc...). The pts do not form a connected graph at 0.7m.
#-----------------------------------------------------------------------------------------
#--linkage map
plot.links(mite.xy, thresh = 0.7)

#--dnearneigh identifies neighbors of region points by Euclidean distance between
#  lower (0) and upper bounds (0.7)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, 0.7)
summary(nb1)

#-----------------------------------------------------------------------------------------
# Correlograms of water density of substrate
#-----------------------------------------------------------------------------------------
watr.dens <- mite.env[,2]
watr.correlog <- sp.correlogram(nb1, watr.dens, order = 14, method = 'I', zero.policy = T)
print(watr.correlog, p.adjust.methods = 'holm')
plot(watr.correlog)

#--This graph shows how correlated pairs of observations are as you increase the distance
# (or lag) between them
