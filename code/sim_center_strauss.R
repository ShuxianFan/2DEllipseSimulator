rm(list = ls())
library(spatstat)
X = rmh(model = list(cif = c('hardcore', 'strauss'),
                 par = list(list(beta = 10, hc = 0.05),
                 list(beta = 1, gamma = 0.5, r = 0.1)),
                 w = owin(xrange = c(0,7), yrange = c(0,0.3))))

# fit <- ppm(X, ~ 1, Hybrid(Hardcore(0.03), Strauss(0.07)))

point.coords = data.frame(x = 1000*X$x, y = 1000*X$y)
source(file = "/Users/shuxianfan/Desktop/2DEllipseSimulator/code/sim_ellipse.R")
elldat = getelldat(point.coords)
plotallEllipse(elldat)+geom_hline(yintercept = c(0,300))

