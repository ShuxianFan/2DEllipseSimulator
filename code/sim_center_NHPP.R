rm(list=ls())

#################
library(tidyverse)
library(mvtnorm)
library(msm)
library(MCMCpack)
library(mcmc)
library(MASS)
library(sp)
library(rgeos)
library(raster)
library(geostatsp)
library(mapmisc)
library(ggplot2)
library(matrixcalc)
##################
# spatial raster
# dimensions  : 15, 350, 5250  (nrow, ncol, ncell)
# resolution  : 20, 20  (x, y)
# extent      : 0, 7000, 0, 300  (xmin, xmax, ymin, ymax)
myRaster = squareRaster(raster(extent(0,7000,0,300), ncols=350,nrows=15))

n = ncell(myRaster) # get number of cells in the raster
raster.coords = xyFromCell(myRaster, 1:n) # get coordinates of the center of raster cells 
distmat = as.matrix(dist(raster.coords)) # euclidean distance matrix 

# generate values RF U(s) using:
# (1)exponential correlation
# (2) gaussian correlation
# (3) Matern correlation - model=c(variance=3,range=0.2, shape=1)

# (1) Exponential variance matrix
phi = 4 
sigma = 1
varmat = sigma^2*exp(-distmat/phi)


# (2) Matern variance matrix
# phi = 0.15
# sigma = 3
# kappa = 1
# varmat = (sigma^2/gamma(kappa))*2^(kappa - 1)*(sqrt(8*kappa)*distmat/phi)^kappa*besselK(sqrt(8*kappa)*distmat/phi, kappa)
# diag(varmat) = sigma^2


# Random resource covariate
thechol = t(chol(varmat)) # get choleski decomposition
U = thechol %*% rnorm(n)
# W = exp(U)
values(myRaster) <- U
plot(myRaster)

# create a new RasterLayer with a higher resolution (smaller cells)
myRaster = disaggregate(myRaster, fact = 2, method = 'bilinear' )
plot(myRaster)
values(myRaster) <- scale(values(myRaster)) # center and scale the covariate

######################
# Simulate event location

beta = matrix(c(-1,1),2,1)  # set parameters
X = cbind(1, values(myRaster)) # design matrix with covariates
ncolX = ncol(X)
lambda = exp(X%*%beta)

# simulate by rejection
n = 1000 # number of proposed points
sample.idx = sample(1:ncell(myRaster), n, replace = T)
keeper = rbinom(n, 1, prob = lambda[sample.idx]/max(lambda))
point.coords = xyFromCell(myRaster, sample.idx[keeper==1])
n.point = nrow(point.coords)
plot(myRaster,xlim = c(0,7000),ylim = c(0, 300))
points(point.coords, col = 1, pch = 18, cex = 0.5)

point.coords <- as.data.frame(point.coords)
colnames(point.coords) <- c('x','y')
source(file = "/Users/shuxianfan/Desktop/2DEllipseSimulator/code/sim_ellipse.R")
elldat = getelldat(point.coords)
plotallEllipse(elldat)+geom_hline(yintercept = c(0,300))


