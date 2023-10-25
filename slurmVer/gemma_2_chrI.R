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

# Genotype Matrix - X
readr::read_csv(system.file("extdata", "mouse100.geno.txt", package = "gemma2"), col_names = FALSE) -> genoT
# Col 1: SNP, Major allele, minor allele, then 1/0 coding for each ind. 
read_csv("./gemma2_test/traits.csv", col_names = FALSE) -> geno
genoI <- geno %>% 
  filter(grepl("I:",X1)) %>% 
  filter(!grepl("II:|III:", X1))
genoI[genoI == 2] <- 1

# Phenotype Matrix - Y
readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE) -> phenoT

read_tsv("./gemma2_test/full_traitData.tsv", col_names = TRUE) -> pheno

hist(pheno$slope)
qqnorm(pheno$slope)
qqline(pheno$slope, col = "red")

# Kinship Matrix - K
readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100] -> kinship
# I would get this matrix from the gemma functionality

read_tsv("./gemma2_test/gemma/grm/gemmeGRM.I.cXX.txt", col_names = FALSE) -> kinshipI

g <- genoI[, - c(1:3)]
t(as.matrix(g)) -> gm # first 5 variants
as.matrix(pheno[ ,c(2,3)]) -> phe23 
as.matrix(pheno[ ,2]) -> phe2 

e_outI <- eigen2(kinshipI)
center_kinship(as.matrix(kinshipI)) -> kinshipI_centered
e_outI <- eigen2(kinshipI_centered)

kinship_vec <- as.matrix(e_outI$vectors)
t(gm) %*% kinship_vec -> XU # note this is X 5 U, NOT X 1 U

MphEM(eval = e_outI$values, 
      X = XU, 
      Y = t(phe2) %*% e_outI$vectors, 
      V_g = matrix(c(0.00754107), nrow = 1), # these starting values came from running GEMMA for one iteration
      V_e = matrix(c(0.00143981), nrow = 1) # need to change the variance components accordingly, properly scale them...test with a 
) -> ground_control

library(Matrix)
MphCalcP <- function(eval, eps_eval, X_vec, W, Y, V_g, V_e) {
  n_size <- length(eval)
  c_size <- nrow(W)
  d_size <- nrow(V_g)
  dc_size <- d_size * c_size
  
  WHix <- matrix(0, nrow = dc_size, ncol = d_size) # dc x d matrix
  QiWHix <- matrix(0, nrow = dc_size, ncol = d_size) # dc x d matrix
  xPx <- matrix(0, nrow = d_size, ncol = d_size) # d x d matrix
  xPy <- matrix(0, nrow = d_size, ncol = 1) # d vector
  WHiy <- matrix(0, nrow = dc_size, ncol = 1) # dc vector
  
  # eigen_proc
  ep_out <- eigen_proc_mod(V_g, V_e)
  ep_out[[1]] -> logdet_Ve # just a quantity
  ep_out[[2]] -> UltVeh # d x d matrix 
  ep_out[[3]] -> UltVehi # d x d matrix
  ep_out[[4]] -> D_l # d vector
  ep_out[[5]] -> V_e_temp # d x d matrix
  # calc_qi
  cq_out <- calc_qi_mod(eval, eps_eval, D_l, W, V_e_temp)
  cq_out[[1]] -> Qi
  cq_out[[2]] -> lndetQ
  
  # Calculate UltVehiY
  UltVehiY <- UltVehi %*% Y
  
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
        dw <- W[j,k]
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
  
  return(p_value)
}

# We'll move on with getting p-values once I know that this officially works
