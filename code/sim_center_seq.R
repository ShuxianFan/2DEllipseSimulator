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
phi = 6
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
# myRaster = disaggregate(myRaster, fact = 2, method = 'bilinear' )
# plot(myRaster)
# values(myRaster) <- scale(values(myRaster)) # center and scale the covariate

######################
# Simulate event location
beta = matrix(c(-1,1),2,1)  # set parameters
X = cbind(1, values(myRaster)) # design matrix with covariates
ncolX = ncol(X)
lambda = exp(X%*%beta)
lambda.raster = squareRaster(raster(extent(0,7000,0,300), ncols=350,nrows=15))
values(lambda.raster) = lambda
######################
# simulate by sequentially updating the lambda 
# simulate by rejection
#  43.33408 417.85153
plot(lambda.raster)

findneighbours = function(sample.idx1, radius){
  ncell = sample.idx1$cell
  distvec= distmat[ncell,]
  neighbours = which(distvec<=radius)
  return(neighbours)
}

updatelambda = function(sample.idx,lambda, radius){
  neighbours = findneighbours(sample.idx1 = sample.idx, radius = radius)
  lambda[neighbours]<-0
  return(lambda)
}


n = 20
i = 1
point.coords = matrix(nrow = 0, ncol = 2)
while (i<=n) {
  sample.idx = as.data.frame(sampleRandom(lambda.raster, 1, cell = TRUE, rowcol = TRUE))
  keeper = rbinom(1, 1, prob = lambda[as.integer(sample.idx$cell)]/max(lambda))
  if(keeper == 1){
    radius = runif(1, min=50, max = 300)
    lambda = updatelambda(sample.idx,lambda,radius)
    values(lambda.raster) <- lambda
    if(i%%2==0){plot(lambda.raster)}
    point.coords = rbind(point.coords,xyFromCell(lambda.raster, sample.idx$cell))
    i = i+1
  }
}

plot(lambda.raster,xlim = c(0,7000),ylim = c(0, 300))
points(point.coords, col = 1, pch = 18, cex = 0.5)

point.coords <- as.data.frame(point.coords)
colnames(point.coords) <- c('x','y')
source(file = "/Users/shuxianfan/Desktop/2DEllipseSimulator/code/sim_ellipse.R")
elldat = getelldat(point.coords)
plotallEllipse(elldat)+geom_hline(yintercept = c(0,300))






