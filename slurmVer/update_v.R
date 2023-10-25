library(readxl)
library(plyr)
library(tidyverse)
library(broom)
#library(ggpubr)

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
load(file = "./gemma2_test/epsilon.Rdata")

# update_v
n_size <- length(eval)
d_size <- nrow(U_l)
V_g <- matrix(0, nrow = d_size, ncol = d_size)
V_e <- V_g
for (k in 1:n_size){
  delta <- eval[k]
  epsilon <- eps_eval[k] / V_e_temp
  #if (delta != 0){
  U_col <- U_hat[, k]
  E_col <- E_hat[, k] # THIS WAS THE PROBLEM 
  V_g <- V_g + U_col %*% t(U_col) / delta
  V_e <- V_e + E_col %*% t(E_col) / epsilon
  #}
}

V_g <- V_g + epsilon_uu
V_e <- V_e + epsilon_ee
V_g <- V_g / n_size
V_e <- V_e / n_size

save(V_g, V_e, file = "./gemma2_test/variances.Rdata")
