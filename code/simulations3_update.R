setwd("/Users/mac/Documents/single_cell_cancer/single_cell_cancer")
library(Matrix)
library(flashier)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

## load in functions
source("code/util.R")

### load in simulated loadings and factors
L <- readRDS("data/model.rds")$L
F <- readRDS("data/model.rds")$F

### set seed
set.seed(100)

### load in simulated loadings and factors
L <- readRDS("data/model.rds")$L
F <- readRDS("data/model.rds")$F
F[,1] <- F[,1] - 10


### sigma2=1
### simulate data with i.i.d. gaussian errors
X <- L %*% t(F) + matrix(rnorm(800*2000, 0, 1), nrow=800, ncol=2000)

### calculate the rescaled XX'
XXt <- X %*% t(X)/ncol(X)

### fit flash with point exponential prior to XXt 
fit3 <- fit.ebcovmf.snn(dat=XXt, Kmax=11)
fit3$sampler <- NULL
fit3$flash.fit <- NULL
saveRDS(fit3, "output/simulations3_update_XXt_fit3_sigma_v1.Rds")


### sigma2=0.5
### simulate data with i.i.d. gaussian errors
X <- L %*% t(F) + matrix(rnorm(800*2000, 0, 0.7), nrow=800, ncol=2000)

### calculate the rescaled XX'
XXt <- X %*% t(X)/ncol(X)

### fit flash with point exponential prior to XXt 
fit3 <- fit.ebcovmf.snn(dat=XXt, Kmax=11)
fit3$sampler <- NULL
fit3$flash.fit <- NULL
saveRDS(fit3, "output/simulations3_update_XXt_fit3_sigma_v2.Rds")