setwd("/scratch/midway2/yushaliu/single_cell_cancer/simulations")
library(Matrix)
library(flashier)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

## load in functions
source("util.R")

### load in simulated loadings and factors
L <- readRDS("model.rds")$L
F <- readRDS("model.rds")$F
F[,1] <- F[,1] - 10

### scaling factors for L and F 
sf <- sqrt(nrow(F)/colSums(F^2))
init.F <- t(t(F)*sf)
init.L <- t(t(L)/sf)



############################################################## sigma2=1 ##############################################################
### run 10 replicates
for(rep in 1:10){
  ### simulate data
  set.seed(100*rep)
  X <- L %*% t(F) + matrix(rnorm(800*2000, 0, 1), nrow=800, ncol=2000)
  
  ### calculate the rescaled XX'
  XXt <- X %*% t(X)/ncol(X)
  
  ### fit flash without considering the diagonal component for now
  fit1 <- flash.init(XXt, var.type = 0) %>% flash.init.factors(
    EF=list(init.L, init.L), prior.family=as.prior(ebnm::ebnm_point_exponential, sign = 1)) 
  ### fit flash again with the diagonal component
  fit1 <- fit.ebcovmf(dat=XXt, fl=fit1, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1))$fl
  saveRDS(fit1, paste0("output/simulations3_more_XXt_fit1_sigma_1_rep", rep, ".rds"))
  
  ### fit flash with point exponential prior to XXt while adding all factors then backfitting
  fit2 <- fit.ebcovmf.snn.v1(dat=XXt, Kmax=11)
  saveRDS(fit2, paste0("output/simulations3_more_XXt_fit2_sigma_1_rep", rep, ".rds"))
  
  ### fit flash with point exponential prior to XXt while adding and backfitting factors iteratively
  fit3 <- fit.ebcovmf.snn(dat=XXt, Kmax=11)
  saveRDS(fit3, paste0("output/simulations3_more_XXt_fit3_sigma_1_rep", rep, ".rds"))
}



print(sessionInfo())