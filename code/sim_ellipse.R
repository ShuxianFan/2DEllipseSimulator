library(tidyverse)
library(mvtnorm)
library(msm)
library(MASS)
library(sp)
library(rgeos)
library(raster)
library(geostatsp)
library(mapmisc)
library(ggplot2)
library(matrixcalc)
######################
# Simulate ellipse shape

# point.coords<-as.data.frame(point.coords)
# colnames(point.coords) = c("x", "y")
# ggplot(point.coords) + geom_point(aes(x,y)) + xlim(0,7000)+ylim(0,300)
# point.dist = as.matrix(dist(point.coords))

# Helper function for generating a random size for the ellipse 
# OUTPUT: a vector of varx, vary, cov for the covariance of ellipse
surf0_dat = read.csv("/Users/shuxianfan/Desktop/2DEllipseSimulator/data/surf0_dat.csv")

randomSize = function(){
  quantiles_varx = quantile(surf0_dat$var_x, probs = c(0, 0.1, 0.9, 1))
  quantiles_vary = quantile(surf0_dat$var_y, probs = c(0, 0.1, 0.9, 1))
  quantiles_cov = quantile(surf0_dat$cov, probs = c(0, 0.1, 0.9, 1))
  
  repeat {
    varx = sample(x = c(runif(1,quantiles_varx[1],quantiles_varx[2]), 
                        runif(1,quantiles_varx[2],quantiles_varx[3]),
                        runif(1,quantiles_varx[3],quantiles_varx[4])),
                  size = 1, prob = c(0.1, 0.8, 0.1))
    
    vary = sample(x = c(runif(1,quantiles_vary[1],quantiles_vary[2]), 
                        runif(1,quantiles_vary[2],quantiles_vary[3]),
                        runif(1,quantiles_vary[3],quantiles_vary[4])),
                  size = 1, prob = c(0.1, 0.8, 0.1))
    
    cov = sample(x = c(runif(1,quantiles_cov[1],quantiles_cov[2]), 
                       runif(1,quantiles_cov[2],quantiles_cov[3]),
                       runif(1,quantiles_cov[3],quantiles_cov[4])),
                 size = 1, prob = c(0.1, 0.8, 0.1))
    # check is positive definite
    if(is.positive.definite(matrix(c(varx,cov,cov,vary),2,2))) break
  }
  return(data.frame(var_x = varx,var_y = vary,cov = cov))
}




# Helper function for deciding if any point is inside a ellipse
# INPUT: a point vector or a matrix of points
# OUTPUT: TRUE if one of the points is inside the ellipse/ FALSE if all points are outside ellipse
is.inside = function(points, ctr,ellCov){
  indicator = FALSE
  if(is.null(nrow(points))){n = length(points)}
  else{n = nrow(points)}
  #  browser()
  for (i in 1:n) {
    if(is.null(nrow(points))){
      x = points[1]
      y = points[2]
    }
    else{
      x = points[i,][1]
      y = points[i,][2]
    }
    A = matrix(c(as.numeric(ellCov['var_x']), as.numeric(ellCov['cov']),
                 as.numeric(ellCov['cov']),
                 as.numeric(ellCov['var_y'])),byrow = T,2,2)
    eigenA = eigen(A)
    major = 2*sqrt(max(eigenA$values)) # major axis
    minor = 2*sqrt(min(eigenA$values)) # minor axis
    rotation = diag(eigenA$vectors)[1] # in radian 
    cos_angle = cos((180-rotation)*pi/180)
    sin_angle = sin((180-rotation)*pi/180)
    xc = x-ctr[1]
    yc = y-ctr[2]
    xct = xc*cos_angle - yc*sin_angle
    yct = xc*sin_angle + yc*cos_angle
    
    rad_cc = (xct^2/(minor/2)^2) + (yct^2/(major/2)^2)
    if(rad_cc<=1){
      indicator <- TRUE
      break
    }
    else{next}
  }
  return(indicator)
}

# Helper function for deciding if two ellipse intersect
# INPUT: centers and covariance of two ellipses
# OUTPUT: TRUE if they intersect/ FALSE if they do not intersect
is.intersect = function(ctr1, ellCov1, ctr2, ellCov2){
  # covariance matrix for ellipse 1 and 2
  #  browser()
  A1 = matrix(c(as.numeric(ellCov1['var_x']), 
                as.numeric(ellCov1['cov']),
                as.numeric(ellCov1['cov']),
                as.numeric(ellCov1['var_y'])),2,2)
  A2 = matrix(c(as.numeric(ellCov2['var_x']), 
                as.numeric(ellCov2['cov']),
                as.numeric(ellCov2['cov']),
                as.numeric(ellCov2['var_y'])),2,2)
  # cholesky decomposition
  RR1 = chol(A1)
  RR2 = chol(A2)
  # angles for ellipse
  angles <- seq(0, 2*pi, length.out=200)          
  ell1    <- 1 * cbind(cos(angles), sin(angles)) %*% RR1 
  ell2   <- 1 * cbind(cos(angles), sin(angles)) %*% RR2  # ellipse scaled with factor 1
  # sweep points to ctr 
  ellCtr1 <- sweep(ell1, 2, ctr1, "+")  
  ellCtr2 <- sweep(ell2, 2, ctr2, "+")  
  # check if ellipse 1 is in ellipse 2 or if ellipse 2 is in ellipse 2
  
  return(is.inside(ellCtr1, ctr2,ellCov2) || is.inside(ellCtr2, ctr1,ellCov1))
}

# Sequentially generating the ellipse for the point centers
# INPUT: point.coords
# OUTPUT: a sequence of ellipses(varx,vary,cov) that not overlapping for each point

# sort the point.coords by x
point.coords<-point.coords%>% arrange(x)

# sanity check
min(dist(point.coords))>50


# create elldat

getelldat = function(point.coords){
  elldat = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(elldat) = c("x","y","var_x","var_y","cov")
  elldat[1,]<- as.vector(c(point.coords[1,], randomSize())) # add first ellipse in 
  i = 2
  while(i<=nrow(point.coords)){
    ctr = c(point.coords$x[i], point.coords$y[i])
    ellCov = randomSize()
    # check if all non-overlap
    #  browser()
    if(all(apply(elldat, 1, function(ell){!is.intersect(ctr1 = ctr, 
                                                        ellCov1 = ellCov, 
                                                        ctr2 = ell[c('x','y')],
                                                        ellCov2 = ell[c('var_x','var_y','cov')])}))){
      elldat[i,] <- as.numeric(c(ctr,ellCov))
      i = i+1
    }
  }
  return(elldat)
}


# elldat = getelldat(point.coords = point.coords)

# Helper function for plotting ellipse 
getell = function(x,y,varx,vary,cov){
  ctr = c(x,y)
  A = matrix(c(varx,cov,cov, vary),2,2)
  RR = chol(A)
  angles <- seq(0, 2*pi, length.out=200)
  ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR
  ellCtr <- sweep(ell, 2, ctr, "+")
  return(data.frame(x = ellCtr[,1], y = ellCtr[,2]))
}

# Helper function for plotting all ellipse on the dimensional surface
plotallEllipse<-function(elldat){
  p <- ggplot(data =getell(x = elldat$x[1], y = elldat$y[1], varx = elldat$var_x[1], vary = elldat$var_y[1], cov = elldat$cov[1]), aes(ymin=0, ymax=300, xmin=0, xmax=7000))+ geom_path(aes(x,y))
  
  for(i in 2:nrow(elldat)){ 
    p <- p+geom_path(aes(x,y), 
                     data = getell(x = elldat$x[i], y = elldat$y[i], varx = elldat$var_x[i], vary = elldat$var_y[i], cov = elldat$cov[i]))
  }
  return(p)
}

# plotallEllipse(elldat)



