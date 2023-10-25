library(readxl)
library(plyr)
library(tidyverse)
library(broom)
library(ggpubr)

#install.packages("factoextra")
library(factoextra)
library(gemma2)

#setwd("/Users/jamesonblount/gemma2_test/")
setwd("/hpc/group/baughlab/jb621/GEMMA_2_testing")
load(file = "./gemma2_test/starter.Rdata")
load("./gemma2_test/randGenoI.eigen.Rdata")
load(file="./gemma2_test/calcQ_f.Rdata")
load(file = "./gemma2_test/XHiY_f.Rdata")

# calc_omega(eval, eps_eval, D_l)
n_size <- length(eval)
d_size <- length(D_l)
OmegaU <- matrix(nrow = d_size, ncol = n_size)
OmegaE <- OmegaU
for (k in 1:n_size){
  delta <- eval[k]
  epsilon <- eps_eval[k] / V_e_temp
  for (i in 1:d_size){
    dl <- D_l[i]
    d_u <- dl / (dl * delta + epsilon)
    d_e <- delta * d_u
    OmegaU[i, k] <- d_u
    OmegaE[i, k] <- d_e
  }
}
save(OmegaU, OmegaE, file = "./gemma2_test/Omega.Rdata")

# UpdateRL_b(xHiy, Qi, d_size = nrow(Y))
d_size <- nrow(Y)
nrow(Qi) -> dc_size
c_size <- dc_size / d_size
b <- Qi %*% xHiy
UltVehiB <- matrix(nrow = d_size, ncol = c_size)
for (i in 1:c_size){
  b_subcol <- b[(1 + (i - 1) * d_size):(i * d_size)]
  b_subcol -> UltVehiB[, i] # could use as.matrix here
}

UltVehiBX <- UltVehiB %*% X

save(b, UltVehiB, file = "./gemma2_test/UpdateRL_b.Rdata")

# update_u
UltVehiU <- UltVehiY
UltVehiU <- UltVehiU - UltVehiBX
UltVehiU <- UltVehiU * OmegaE

# update_e
UltVehiE <- UltVehiY - UltVehiBX - UltVehiU

U_hat <- t(UltVeh) %*% UltVehiU
E_hat <- t(UltVeh) %*% UltVehiE
B <- t(UltVeh) %*% UltVehiB

save(UltVehiU, UltVehiE, U_hat, E_hat, B, file = "./gemma2_test/Update_UE.Rdata")
