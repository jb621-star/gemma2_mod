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
load(file = "./gemma2_test/epsilon.Rdata")
load(file = "./gemma2_test/variances.Rdata")

# Before going any further, checking to see how these quantities differ
MphEM(eval = e_outI$values, 
      X = XU, 
      Y = t(phe2) %*% e_outI$vectors, 
      V_g = matrix(c(0.00754107), nrow = 1), # these starting values came from running GEMMA for one iteration
      V_e = matrix(c(0.00143981), nrow = 1) # need to change the variance components accordingly, properly scale them...test with a 
) -> ground_control
  

# Calculate UltVehiY
UltVehiY <- UltVehi %*% Y

# As input for this function, there is
# X_row.vector <- gsl_vector_view X_row = gsl_matrix_row(X, c_size); a vector of the transformed marker info
# X_sub.matrix <- gsl_matrix_transpose_memcpy(&X_sub.matrix, UtW); copied UtW
# What is UtW? 
# X_vec is a given row vector from X
# W is a "sub-matrix of X", that being the transposed copy of UtW, which is eigenvectors times covariates
W <- matrix(1, nrow = length(Y))
read_tsv("./gemma2_test/gemma/grm/gemmeGRM.I.cXX.txt", col_names = FALSE) -> kinshipI
e_outI <- eigen2(kinshipI)
center_kinship(as.matrix(kinshipI)) -> kinshipI_centered
e_outI <- eigen2(kinshipI_centered)

t(W) %*% e_outI$vectors -> UtW

n_size <- length(eval)
c_size <- nrow(UtW) # for this case should be 100, so we're getting p-values one SNP at a time
d_size <- nrow(V_g)
dc_size <- d_size * c_size

WHix <- matrix(0, nrow = dc_size, ncol = d_size) # dc x d matrix
QiWHix <- matrix(0, nrow = dc_size, ncol = d_size) # dc x d matrix
xPx <- matrix(0, nrow = d_size, ncol = d_size) # d x d matrix
xPy <- matrix(0, nrow = d_size, ncol = 1) # d vector
WHiy <- matrix(0, nrow = dc_size, ncol = 1) # dc vector
p_mat <- matrix(0,nrow = nrow(X), ncol = 1) # vector with one row per SNP

for (h in 1:nrow(X)) {
  X_vec <- X[h,]


# Calculate WHix, WHiy, xHiy, xHix
for (i in 1:d_size){
  dl <- D_l[i]
  
  d1 = 0
  d2 = 0
  
  for (k in 1:n_size){
    dx <- X_vec[k]
    dy <- UltVehiY[i, k]
    delta <- eval[k]
    epsilon <- eps_eval[k] / V_e_temp
    
    d1 <- d1 + dx * dy / (delta * dl + epsilon)
    d2 <- d2 + dx * dx / (delta * dl + epsilon)
  }
  xPy[i] <- d1
  xPx[i,i] <- d2
  
  for (j in 1:c_size) {
    d1 = 0
    d2 = 0
    for (k in 1:n_size) {
      delta <- eval[k]
      epsilon <- eps_eval[k] / V_e_temp
      dx <- X_vec[k]
      dw <- UtW[j,k]
      dy <- UltVehiY[i,k]
      
      d1 <- d1 + dx * dw / (delta * dl + epsilon);
      d2 <- d2 + dy * dw / (delta * dl + epsilon);
    }
    WHix[i, i] <- d1
    WHiy[i] <- d2
  }
}
QiWHix <- Qi %*% WHix
xPx <- WHix %*% QiWHix
xPy <- QiWHix %*% WHiy

# Calculate V(beta) and beta
# LU decomposition
lu_decomp <- lu(xPx)

# Solve linear system
D_l <- solve(lu_decomp@x, xPy)

# Compute matrix inverse
Vbeta <- solve(xPx)

# Need to multiply UltVehi on both sides or one side
# Matrix multiplication
beta <- t(UltVeh) %*% D_l
xPx <- Vbeta %*% UltVeh
Vbeta <- t(UltVeh) %*% xPx

# Dot product to calculate test statistic and pvalue
d <- sum(D_l * xPy)

#p-value from test statistic
p_value <- 1 - pchisq(d, d_size)
p_mat[h] <- p_value
}
save(p_mat, file = "./gemma2_test/pvals.Rdata")
