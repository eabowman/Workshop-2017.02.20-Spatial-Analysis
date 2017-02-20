## Script created by Liz Bowman for use in a workshop on how to examine
## data for spatial paterns and to test for correlations between communities
## and environmental variables
## eabowman@email.arizona.edu
## created on Feb. 20, 2017

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

#--if you are interested in what all of the data is then check out the help page
#?mite

#=========================================================================================
# Spatial correlogram based on Moran's I
#=========================================================================================

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
#--isolate the column with water density data
watr.dens <- mite.env[,2]
#--make a spatial correlogram
watr.correlog <- sp.correlogram(nb1, watr.dens, order = 14, method = 'I', zero.policy = T)
#--print the p-values for each distance class
print(watr.correlog, p.adjust.methods = 'holm')
#--plot spatial pattern
plot(watr.correlog)

#--This graph shows how correlated pairs of observations within each distance class are
# as distance increases (or lag) between them
