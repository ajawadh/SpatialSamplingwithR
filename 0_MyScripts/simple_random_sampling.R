library(remotes)
library(sswr)
library(ggplot2)
library(dplyr)
library(QuantileNPCI)

# SRSWOR from discretized grids with points at the center, each square is 25 m in length
n <- 40
N <- nrow(grdVoorst) # grdVoorst is organized in a discritiyed grid in a tibble
# where each square center is saved. this is the sampling frame.
#set.seed(314)
#units <- sample(N, size = n, replace = FALSE)
#mysample <- grdVoorst[units, ] 

# print(mysample)
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

# Check sampled indices
print(units)

mysample <- grdVoorst[units, ]

cellsize <- 25
mysample$s1 <- jitter(mysample$s1, amount = cellsize / 2)
mysample$s2 <- jitter(mysample$s2, amount = cellsize / 2)

# calculating the mean of sample
mz <- mean(mysample$z)

# multiplying this by the area alone is not useful (different units)
# calculated the volumen of soil in the grids for the first 0-30 cm
# each grid is multiplied by the number of grids N * size of grid * thikness of soil
# total SOC is then estimated by volume * bulk density of soil * SOC mean
# this is multiplied by 10-6 to obtain total mass of SOC in Mg (1000 kg)

vol_soil <- N * 25^2 * 0.3
bd <- 1500
tz <- vol_soil * bd * mz * 10^-6
print(tz)

# plotting the cdf
dev.new()
p <- ggplot(mysample, mapping = aes(z)) +
     stat_ecdf(geom = "step") +
     scale_x_continuous(name = "SOM") +
     scale_y_continuous(name = "F")

print(p)

quantile(mysample$z, probs = c(0.25, 0.5, 0.75), type = 4) %>%
  round(1)

res <- quantCI(mysample$z, q = 0.5, alpha = 0.05, method = "exact")
print(res)

se_mz <- sqrt(var(mysample$z) / n)
print(se_mz)

se_tz <- se_mz * vol_soil * bd * 10^-6
print(se_tz)


library(survey)
mysample$pi <- n / N
design_si <- svydesign(id = ~ 1, probs = ~ pi, data = mysample)
print(svymean(~ z, design = design_si))

mysample$N <- N
design_si_fp <- svydesign(id = ~ 1, probs = ~ pi, fpc = ~ N, data = mysample)
print(svymean(~ z, design_si_fp))

print(svyquantile(~ z, design_si, quantile = c(0.5, 0.9)))

alpha <- 0.05
margin <- qt(1 - alpha / 2, n - 1, lower.tail = TRUE) * se_mz
lower <- mz - margin
upper <- mz + margin

confint(svymean(~ z, design_si), df = degf(design_si), level = 0.95)