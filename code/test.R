#########################################

# The idea basically is: 
# the major and minor half-diameters are the two eigen values 
# and you rotate the ellipse by the amount of angle between
# the first eigen vector and the x-axis

ctr    <- c(surf0_dat$x[1], surf0_dat$y[1])                              # data centroid -> colMeans(dataMatrix)
A      <- matrix(c(surf0_dat$var_x[1],surf0_dat$cov[1],surf0_dat$cov[1],surf0_dat$var_y[1]), nrow=2) # covariance matrix -> cov(dataMatrix)
RR     <- chol(A)                               # Cholesky decomposition
angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
plot(ellCtr, type="l", lwd=2, asp=1)            # plot ellipse
points(ctr[1], ctr[2], pch=4, lwd=2)            # plot data centroid

hist(surf0_dat$var_x)
hist(surf0_dat$var_y)
hist(surf0_dat$cov)

boxplot(surf0_dat$var_x)
boxplot(surf0_dat$var_y)
boxplot(surf0_dat$cov)


randomSize = function(surf0_dat){
  quantiles_varx = quantile(surf0_dat$var_x, probs = c(0, 0.1, 0.9, 1))
  quantiles_vary = quantile(surf0_dat$var_y, probs = c(0, 0.1, 0.9, 1))
  quantiles_cov = quantile(surf0_dat$cov, probs = c(0, 0.1, 0.9, 1))
  
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
  return(c(varx,vary,cov))
}


randomSize(surf0_dat)

# Input: center, varx, vary, cov 
# output: ellipse 
# The idea basically is: 
# the major and minor half-diameters are the two eigen values 
# and you rotate the ellipse by the amount of angle between
# the first eigen vector and the x-axis

ctr    <- c(surf0_dat$x[1], surf0_dat$y[1])                              # data centroid -> colMeans(dataMatrix)
A      <- matrix(c(surf0_dat$var_x[1],surf0_dat$cov[1],surf0_dat$cov[1],surf0_dat$var_y[1]), nrow=2) # covariance matrix -> cov(dataMatrix)
RR     <- chol(A)                               # Cholesky decomposition
angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
plot(ellCtr, type="l", lwd=2, asp=1)            # plot ellipse
points(ctr[1], ctr[2], pch=4, lwd=2)            # plot data centroid

eigenA = eigen(A)

major = 2*sqrt(max(eigenA$values)) # major axis
minor = 2*sqrt(min(eigenA$values)) # minor axis
rotation = diag(eigenA$vectors)[1]*180/pi



size = randomSize(surf0_dat)
ctr    <- c(surf0_dat$x[1], surf0_dat$y[1])                              # data centroid -> colMeans(dataMatrix)
A      <- matrix(c(size[1],size[3],size[3],size[2]), nrow=2) # covariance matrix -> cov(dataMatrix)
RR     <- chol(A)                               # Cholesky decomposition
angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
plot(ellCtr, type="l", lwd=2, asp=1)            # plot ellipse
points(ctr[1], ctr[2], pch=4, lwd=2)    

eigenA = eigen(A)

major = 2*sqrt(max(eigenA$values)) # major axis
minor = 2*sqrt(min(eigenA$values)) # minor axis
rotation = diag(eigenA$vectors)[1] # in radian 

# radians = degrees * (pi/180)
# degrees = radians * (180/pi)
cos_angle = cos((180-rotation)*pi/180)
sin_angle = sin((180-rotation)*pi/180)

# x y is the point axis

xc = x-ctr[1]
yc = y-ctr[2]

xct = xc*cos_angle - yc*sin_angle
yct = xc*sin_angle + yc*cos_angle

rad_cc = (xct^2/(minor/2)^2) + (yct^2/(major/2)^2)

# if rad_cc <=1, point is in ellipse



