library(ggplot2)
library(tidyverse)
common_path = "/Users/shuxianfan/Desktop/2DEllipseSimulator/data/16Mar2016/"
# 16Mar2016
# 16Oct2015
# 21Oct2015
primary_dirs = list.files(common_path)

getBoardList = function(common_path){
  primary_dirs = list.files(common_path)
  for (dir in primary_dirs) {
    assign(paste("b_",dir,sep = ""),
           read.csv(file = paste(common_path,dir,"/enhanced_matching.csv", sep = "")))
  }
  # get list of boards and add id
  board_list = lapply(primary_dirs, function(dir){
    read.csv(file = paste(common_path,dir,"/enhanced_matching.csv", sep = ""),header=T) %>% 
      mutate(id = paste(dir))})
  
  return(board_list)
}


board_list1 = getBoardList(common_path = "/Users/shuxianfan/Desktop/2DEllipseSimulator/data/16Mar2016/")
board_list2 = getBoardList(common_path = "/Users/shuxianfan/Desktop/2DEllipseSimulator/data/21Oct2015/")


############# Board List Reading Ends Here #########################
names(board_list1) <- paste("lumber", 1:length(board_list1), sep = "")
list2env(board_list1 , envir = .GlobalEnv)
surf3 = lumber1[lumber1$surface==3,]
surf2 = lumber1[lumber1$surface==2,]
surf1 = lumber1[lumber1$surface==1,]
surf0 = lumber1[lumber1$surface==0,]
plot(surf0$x, surf0$y)
plot(surf1$x, surf1$y)
plot(surf2$x, surf2$y)
plot(surf3$x, surf3$y)


# surface 0 and 2 are wider surface of a piece of lumber 
library(spatstat)
surf0.pp = ppp(surf0$x,surf0$y,c(min(surf0$x)-100,max(surf0$x)+100),c(min(surf0$y)-100,max(surf0$y)+100))
plot(Gest(surf0.pp))
plot(Fest(surf0.pp))
plot(Kest(surf0.pp))


# Combining all points from surface 0
points2D = data.frame(matrix(ncol = 2, nrow = 0))
for(i in 1:length(board_list1)){
  board1 = board_list1[[i]]
  surf0 = board1[board1$surface==0,]
  points2D = rbind(points2D, data.frame(x = surf0$x, y = surf0$y))
}
for(i in 1:length(board_list2)){
  board1 = board_list1[[i]]
  surf0 = board1[board1$surface==0,]
  points2D = rbind(points2D, data.frame(x = surf0$x, y = surf0$y))
}

plot(points2D$x, points2D$y, pch = 20, cex = 0.25)
library(dbscan)
library(factoextra)
db = dbscan(points2D,115, minPts = 5)
print(db)

fviz_cluster(db, points2D, geom = "point")

# Check the varx and vary and cov distribution
surfaces0 = lapply(board_list1, function(lumber){lumber[lumber$surface==0,]})
save(surfaces0, file="surf0list.RData")

surf0_dat = Reduce(rbind, surfaces0)
hist(surf0_dat$var_x)
hist(surf0_dat$var_y)
hist(surf0_dat$cov)

mean(surf0_dat$var_x)
sd(surf0_dat$var_x)

mean(surf0_dat$var_y)
sd(surf0_dat$var_y)

mean(surf0_dat$cov)
sd(surf0_dat$cov)

dist.list = lapply(surfaces0, function(sf){dist(sf[,1:2])})
min.dists.surf0  = unlist(lapply(dist.list, function(sf){min(sf)}))
hist(min.dists.surf0)
range(min.dists.surf0)

write.csv(surf0_dat, file = "surf0_dat.csv")



