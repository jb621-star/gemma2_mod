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
g <- distinct(g)
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
max_iter <- 10000
eval <- e_outI$values # In the original math, these are Dk, the eigenvalues from the decomposition of the kinship matrix, where k is the strain
eps_eval <- eps_outI # these will be e, where k is the strain
X <- XU
Y <- t(phe2) %*% e_outI$vectors
V_g <- matrix(c(0.00754107), nrow = 1)
V_e <- matrix(c(0.00143981), nrow = 1)
V_e -> V_e_temp

save(eval, eps_eval, X, Y, V_g, V_e, V_e_temp, file = "./gemma2_test/starter.Rdata")
#load(file = "./gemma2_test/starter.Rdata")
load(file = "./gemma2_test/Qmat.Rdata")
load(file="./gemma2_test/randGenoI.eigen.Rdata")
load(file="./gemma2_test/calcQ_f.Rdata")

UltVehiY <- UltVehi %*% Y

# calc_XHiY(eval, eps_eval, D_l, X, UltVehiY)
n_size <- length(eval)
c_size <- nrow(X)
d_size <- length(D_l)
xHiy <- rep(0, d_size * c_size)
for (i in 1:d_size){
  dl <- D_l[i]
  for (j in 1:c_size){
    d <- 0
    for (k in 1:n_size){
      x <- X[j, k]
      y <- UltVehiY[i, k]
      delta <- eval[k]
      epsilon <- eps_eval[k] / V_e_temp# the eigenvalue from eigen decomp of resid var                                  from the kth strain
      d <- d + x * y / (delta * dl + epsilon)
    }
    xHiy[(j - 1) * d_size + i] <- d
  }
}

save(UltVehiY, xHiy, file = "./gemma2_test/XHiY_f.Rdata")
