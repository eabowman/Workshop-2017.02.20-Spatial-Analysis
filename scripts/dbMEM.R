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
#install.packages('adespatial')
#install.packages('ade4)

#--once above packages are installed run the following code to load them
library(vegan); library(adespatial); library(ade4)

#--source the sr.value function script located in the script folder
source('./scripts/sr.value.R')

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

#--normalize community data so that it can be used with euclidean distance
mite.h <- decostand(mite,'hellinger')

#--test for linear trends; if signficiant, then we will need to detrend the data
anova(rda(mite.h,mite.xy))

#--to detrend, we will regress all variables against the x-y coordinates and retain
#--the residuals
#--computation of linearly detrended mite data
mite.h.det <- resid(lm(as.matrix(mite.h) ~ .,data = mite.xy))

#=========================================================================================
# Distance based Moran's Eigenvector Maps
#=========================================================================================

#-----------------------------------------------------------------------------------------
# 1a. construct the matrix of PCNM variables step by step...
#-----------------------------------------------------------------------------------------

#--distance matrix of spatial data
xy.d1 <- dist(mite.xy)
#--search for the truncation threshold: maximum value of the min. spanning tree of the
#--euclidean distance matrix using vegan's function spanntree. Use that distance, or any
#--other distance larger than that, as the truncation distance.
spanning <- spantree(xy.d1) #--finds a min. spanning tree connecting all pts.
dmin <- max(spanning$dist)
dmin

#--Truncate te distance matrix
xy.d1[xy.d1 > dmin] <- 4*dmin
#--PCoA of truncated distance matrix
xy.PCoA <- cmdscale(xy.d1, k=nrow(mite.xy)-1, eig = T)
#--count pos. eigenvalues
nb.ev <- length(which(xy.PCoA$eig > 0.0000001))
mite.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(mite.xy), 1:nb.ev])

#-----------------------------------------------------------------------------------------
# 2. Run the global PCNM analysis on the *detrended* mite data
#-----------------------------------------------------------------------------------------

mite.PCNM.rda <- rda(mite.h.det, mite.PCNM)
anova.cca(mite.PCNM.rda)

#-----------------------------------------------------------------------------------------
# 3. Since the analysis is signficant, compute the adjusted R2 and run a forward selection
# of the PCNM vairables
#-----------------------------------------------------------------------------------------

mite.R2a <- RsquareAdj(mite.PCNM.rda)$adj.r.squared
#--Performa a forward selection by permutation of residuals under reduced model
# forward.sel(Y, X, adjRthresh); Y = response data matrix, X = explanatory data matrix
# stop the forward selection procedure if the R-square 
mite.PCNM.fwd <- forward.sel(mite.h.det, as.matrix(mite.PCNM), adjR2thresh = mite.R2a)

#--number of significant PCNM
nb.sig.PCNM <- nrow(mite.PCNM.fwd)
#--identity of significant PCNMs in increasing order
PCNM.sign <- sort(mite.PCNM.fwd[,2])
#--write the significant PCNMs to a new object
PCNM.red <- mite.PCNM[c(PCNM.sign)]

#-----------------------------------------------------------------------------------------
# 4. New PCNM analysis with 10 significant PCNM variables
# Adjusted R-square after forward selection: R2adj=0.2713
#-----------------------------------------------------------------------------------------

mite.PCNM.rda2 <- rda(mite.h.det ~ ., data=PCNM.red)
mite.fwd.R2a <- RsquareAdj(mite.PCNM.rda2)$adj.r.squared
anova.cca(mite.PCNM.rda2)
axes.test <- anova.cca(mite.PCNM.rda2, by='axis')
#--number of significant axes
nb.ax <- length(which(axes.test[,4] <= 0.05))

#-----------------------------------------------------------------------------------------
# 5. Plot the two significant canonical axes
#-----------------------------------------------------------------------------------------
#--isolate the site scores for axes 1 and 2 (choices) from the ordination (rda)
mite.PCNM.axes <- scores(mite.PCNM.rda2, choices=c(1,2), display = 'lc', scaling = 1)
#--plot the first 2 significant axes
par(mfrow=c(1,2))
sr.value(mite.xy, mite.PCNM.axes[,1])
sr.value(mite.xy, mite.PCNM.axes[,2])

#--close par function above
dev.off()

#-----------------------------------------------------------------------------------------
# Is the above variation realted to any of the environmental variables?
#-----------------------------------------------------------------------------------------

#--Interpreting the spatial varation: regression of the two significant canonical axes on 
# the environmental variables (after normality tests)

#--axis 1
shapiro.test(resid(lm(mite.PCNM.axes[,1] ~ ., data = mite.env)))
mite.PCNM.axis1.env <- lm(mite.PCNM.axes[,1] ~., data = mite.env)
summary(mite.PCNM.axis1.env)

#--axis 2
shapiro.test(resid(lm(mite.PCNM.axes[,2] ~ ., data = mite.env)))
mite.PCNM.axis2.env <- lm(mite.PCNM.axes[,2] ~., data = mite.env)
summary(mite.PCNM.axis2.env)