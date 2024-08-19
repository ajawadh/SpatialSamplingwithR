library(remotes)
library(sswr)
library(ggplot2)
library(dplyr)
library(QuantileNPCI)

N_h <- tapply(grdVoorst$stratum, INDEX = grdVoorst$stratum, FUN = length)
w_h <- N_h / sum(N_h)
n <- 40
print(n_h <- round(n * w_h))

n_h[1] <- n_h[1] - 1

library(sampling)
ord <- unique(grdVoorst$stratum)
print(ord)

set.seed(314)
units <- sampling::strata(
  grdVoorst, stratanames = "stratum", size = n_h[ord], method = "srswr")
mysample <- getdata(grdVoorst, units) %>%
  mutate(s1 = s1 %>% jitter(amount = 25 / 2),
         s2 = s2 %>% jitter(amount = 25 / 2))

#print(units)
#print(mysample)