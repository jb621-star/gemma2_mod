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

read_csv("./gemma2_test/traits.csv", col_names = FALSE) -> geno
genoI <- geno %>% 
  filter(grepl("I:",X1)) %>% 
  filter(!grepl("II:|III:", X1))

read_tsv("./gemma2_test/full_traitData.tsv", col_names = TRUE) -> pheno
read_tsv("./gemma2_test/gemma/grm/gemmeGRM.I.cXX.txt", col_names = FALSE) -> kinshipI
g <- genoI[, - c(1:3)]
#g <- distinct(g)
g <- unique(g)
t(as.matrix(g)) -> gm # first 100 variants
as.matrix(pheno[ ,c(2,3)]) -> phe23 
as.matrix(pheno[ ,2]) -> phe2 

e_outI <- eigen2(kinshipI)
center_kinship(as.matrix(kinshipI)) -> kinshipI_centered
e_outI <- eigen2(kinshipI_centered)

resid_var <- read_tsv("./gemma2_test/Strain_TotalVar.tsv", col_names = FALSE)
eps_outI <- resid_var$X1 # taking the residuals and scaling them by the average variance Ve
kinship_vec <- as.matrix(e_outI$vectors)

t(gm) %*% kinship_vec -> XU # note this is X 5 U, NOT X 1 U
load("./gemma2_test/randGenoI.eigen.Rdata")

max_iter <- 10000
eval <- e_outI$values # In the original math, these are Dk, the eigenvalues from the decomposition of the kinship matrix, where k is the strain
eps_eval <- eps_outI # these will be e, where k is the strain
X <- XU
Y <- t(phe2) %*% e_outI$vectors
V_g <- matrix(c(0.00754107), nrow = 1)
V_e <- matrix(c(0.00143981), nrow = 1)
V_e -> V_e_temp

n_size <- length(eval)
d_size <- length(D_l) # D_l = 5.98, U_l = 1
c_size <- nrow(X)# what is c_size? it's the number of rows in the transposed genotypes matrix
# transposed genotypes matrix is c by n
dc_size <- d_size * c_size
Q <- matrix(0, nrow = dc_size, ncol = dc_size)
# sum x * xt / V_l , V_l = (D_lambda)
# What is this loop actually calculating?
for (i in 1:c_size){
  for (j in 1:c_size){
    for (l in 1:d_size){
      dl <- D_l[l] # the eigen decomp of lambda
      if (j < i){
        d <- Q[(j - 1) * d_size + l, (i - 1) * d_size + l]
      } else {
        d <- 0
        for (k in 1:n_size){
          d1 <- X[i, k] # d1 <- ith SNP, kth strain
          d2 <- X[j, k] # d2 <- jth SNP, kth strain
          delta <- eval[k] # the eigenvalue from the eigen decomp of kinship from                                 the kth strain
          epsilon <- eps_eval[k] / V_e_temp # the eigenvalue from eigen decomp of resid var                                  from the kth strain
          
          d <- d + d1 * d2 / (dl * delta + epsilon) # the denominator is H
        }
      }
      Q[(i -1) * d_size + l, (j - 1) * d_size + l] <- d # numbers are getting big
      # because instead of basically dividing by 1, we're now dividing by a very small number instead of 1, the residual variance is dominating here
    }
  }
}
#save(Q, file = "./gemma2_test/Qmat.Rdata")
load(file = "./gemma2_test/Qmat.Rdata")
#eigen(Q)$values
Qi <- solve(Q, tol = 1e-30)
detQ <- det(Q)
lndetQ <- log(detQ)

save(Qi,detQ,lndetQ, file="./gemma2_test/calcQ_f.Rdata")