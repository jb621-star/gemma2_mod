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
load(file = "./gemma2_test/Omega.Rdata")
load(file = "./gemma2_test/UpdateRL_b.Rdata")
load(file = "./gemma2_test/Update_UE.Rdata")

# calc_epsilon
n_size <- length(eval)
c_size <- nrow(X)
d_size <- length(D_l)
dc_size <- nrow(Qi)
epsilon_ee <- matrix(0, nrow = d_size, ncol = d_size)
epsilon_uu <- epsilon_ee
for (k in 1:n_size){
  OmegaU_col <- OmegaU[, k]
  OmegaE_col <- OmegaE[, k]
  diag(epsilon_uu) <- diag(epsilon_uu) + OmegaU_col
  diag(epsilon_ee) <- diag(epsilon_ee) + OmegaE_col
}
M_u <- matrix(0, nrow = dc_size, ncol = d_size)
M_e <- M_u
for (k in 1:n_size){
  delta <- eval[k]
  epsilon <- eps_eval[k] / V_e_temp
  for (i in 1:d_size){
    dl <- D_l[i]
    for (j in 1:c_size){
      x <- X[j, k]
      d <- x / (dl * delta + epsilon)
      M_e[(j - 1) * d_size + i, i] <- d
      M_u[(j - 1) * d_size + i, i] <- d * dl
    }
  }
  QiM <- Qi %*% M_u
  epsilon_uu <- epsilon_uu + t(M_u) %*% QiM * delta
  QiM <- Qi %*% M_e
  epsilon_ee <- epsilon_ee + t(M_e) %*% QiM * epsilon
} # end loop in k
M <- epsilon_uu %*% UltVeh
epsilon_uu <- t(UltVeh) %*% M
M <- epsilon_ee %*% UltVeh
epsilon_ee <- t(UltVeh) %*% M

save(epsilon_uu, epsilon_ee, file = "./gemma2_test/epsilon.Rdata")