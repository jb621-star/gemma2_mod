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

n_size <- length(eval)
c_size <- nrow(X)
d_size <- nrow(Y)
dc_size <- c_size * d_size
XXt <- X %*% t(X)
logl_const <- -(n_size - c_size) * d_size * log(2 * pi)/2 + d_size * log(det(XXt))/2

load("./gemma2_test/randGenoI.eigen.Rdata")
load(file="./gemma2_test/calcQ_f.Rdata")
load(file = "./gemma2_test/XHiY_f.Rdata")

# Mph CalcLogL(eval, sig_eval, xHiy, D_l, UltVehi, Qi)
n_size <- length(eval)
d_size <- length(D_l) # d is number of phenotypes
dc_size <- nrow(Qi)
logl <- 0
for (k in 1:n_size){
  delta <- eval[k]
  epsilon <- eps_eval[k] / V_e_temp
  for (i in 1:d_size){
    y <- UltVehiY[i, k]
    dl <- D_l[i]
    d <- (dl * delta + epsilon)
    logl <- logl + y^2 / d + log(d)
  }
}
Qiv <- Qi %*% xHiy
d <- t(xHiy) %*% Qiv
stopifnot(length(d) == 1)
logl <- logl - d
MphCalcLogL_edit <- (- 0.5 * logl)

logl_new <- logl_const + MphCalcLogL_edit - 0.5 * n_size * logdet_Ve
logl_new <- logl_new - 0.5 * (lndetQ - c_size * logdet_Ve)
logl_old <- logl_new

save(logl_const, MphCalcLogL_edit, logl_old, file = "./gemma2_test/MphCalcLogL.Rdata")