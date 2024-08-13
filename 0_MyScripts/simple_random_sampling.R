library(remotes)
library(sswr)

# SRSWOR from discretized grids with points at the center, each square is 25 m in length
n <- 40
N <- nrow(grdVoorst) # grdVoorst is organized in a discritiyed grid in a tibble
# where each square center is saved. this is the sampling frame.
set.seed(314)
units <- sample(N, size = n, replace = FALSE)
mysample <- grdVoorst[units, ] 

print(mysample)
# A vector with the centers of the selected cells of the discretisation grid, 
# referred to as discretisation points.
# To void restricting the sampling units to the discretisation points, 
# a simple random sample of points is selected in two stages. First,
# n times a grid cell is selected by SRSWR.
# Second, everytime a grid cell is selected, one point is selected fully
# randomly from that grid cell. This account for the infinite number of points
# in the population. This second step of this selection procedure is implemented
# with function "jitter". it adds random noise to the spatial coordinates of grid centers
# by drawing a continous uniform distribution unif(-c,c) with c half the side length of 
# square grid cells. Compared to the above selection, here we respect that the 
# population is actually infinite.
set.seed(314)
units <- sample(N, size=n, replace = TRUE)
my_sample <- grdVoorst[units, ]
cellsize <- 25
mysample$s1 <- jitter(mysample$s1, amount = cellsize / 2)
mysample$s2 <- jitter(my_sample$s2, amount = cellsize / 2)
print(mysample)
?sswr