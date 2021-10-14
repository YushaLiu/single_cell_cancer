## initialize the snn fit using flash fit without constraints on L
init.snn.LL <- function(dat) {
  LL <- dat$flash.fit$EF[[1]]
  FF <- dat$flash.fit$EF[[2]]
  
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  FF <- cbind(FF, -FF)
  
  to.keep <- (colSums(LL) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  
  return(list(LL, FF))
}


### apply the snn to a data matrix with nonnegative constraints on L
fit.snn <- function(data, Kmax, tol=.Machine$double.eps){
  fit.flashier.pn <- flash.init(data, S = 1/sqrt(nrow(data)), var.type = 2) %>% flash.add.greedy(Kmax = Kmax, 
  prior.family = prior.point.normal(), tol=tol, verbose.lvl = 1)
  
  fit.flashier.pn <- flash.init(data, S = 1/sqrt(nrow(data)), var.type = 2
  ) %>% flash.init.factors(EF = fit.flashier.pn$flash.fit$EF, EF2 = fit.flashier.pn$flash.fit$EF2, prior.family = prior.point.normal()
  ) %>% flash.backfit(verbose.lvl = 1)
  
  snn.init <- init.snn.LL(fit.flashier.pn)
  
  fit.flashier.snn <- flash.init(data, S = 1/sqrt(nrow(data)), var.type = 2
  ) %>% flash.init.factors(EF = snn.init, prior.family = c(prior.nonnegative(), prior.normal.scale.mix())
  ) %>% flash.backfit(verbose.lvl = 1)
  
  kset <- (length(fit.flashier.snn$pve) - rank(fit.flashier.snn$pve) < Kmax) & (fit.flashier.snn$pve > 0)
  
  fit.flashier.snn <- flash.init(data, S = 1/sqrt(nrow(data)), var.type = 2
  ) %>% flash.init.factors(EF = lapply(fit.flashier.snn$flash.fit$EF, function(x) x[, kset]), 
                           EF2 = lapply(fit.flashier.snn$flash.fit$EF2, function(x) x[, kset]),
                           prior.family = c(as.prior(ebnm::ebnm_point_exponential, sign = 1), prior.normal.scale.mix())
  ) %>% flash.backfit(verbose.lvl = 1)
  
  fit.flashier.snn$sampler <- NULL
  fit.flashier.snn$flash.fit <- NULL
  
  return(fit.flashier.snn)
}


## initialize the nonnegative fit to covariance matrix XX' s.t. E[XX'] = LL'+ D based on unconstrained estimate of L
init.snn.cov <- function(fl, kset=1:ncol(fl$flash.fit$EF[[1]])) {
  LL <- fl$flash.fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash.fit$EF[[2]][, kset, drop = FALSE]
  
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash.fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash.fit$EF[[2]][, -kset, drop = FALSE], FF)
  
  to.keep <- (colSums(LL) > .Machine$double.eps) & (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  
  return(list(LL, FF))
}


### apply flash to covariance matrix XX' s.t. E[XX'] = LL'+ D, where D = sigma2*I
fit.ebcovmf <- function(dat, fl, prior=prior.point.laplace(), verbose=1){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf
  
  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    kset <- fl$pve > 0
    fl <- flash.init(dat_minuss2) %>% flash.init.factors(EF = lapply(fl$flash.fit$EF, function(x) x[, kset, drop = FALSE]),
      EF2 = lapply(fl$flash.fit$EF2, function(x) x[, kset, drop = FALSE]), prior.family = prior) %>% flash.backfit(verbose.lvl = verbose)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat=dat, fl=fl, s2=s2))
}


### apply flash to covariance matrix XX' s.t. E[XX'] = LL'+ D, where D = sigma2*I, and only a specified subset of factors are updated
fit.ebcovmf.kset <- function(dat, fl, prior=prior.point.laplace(), kset=1:ncol(fl$flash.fit$EF[[1]])){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf
  
  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    fl <- flash.init(dat_minuss2) %>% flash.init.factors(EF = fl$flash.fit$EF, EF2 = fl$flash.fit$EF2, prior.family = prior) %>% flash.backfit(
      kset = kset, verbose.lvl = 0)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat=dat, fl=fl, s2=s2))
}


### run flash to the covariance matrix by adding factors greedily one by one, before backfitting all existing nonnegative factors
fit.ebcovmf.snn <- function(dat, Kmax){
  # fit unconstrained flash for K=1 without considering the diagonal component 
  fit.cov <- flash.init(dat, var.type = 0) %>% flash.add.greedy(Kmax = 1, prior.family = prior.point.laplace(), verbose.lvl = 1
  ) %>% flash.backfit(verbose.lvl = 1)
  
  # initialize the snn fit based on the unconstrained flash fit 
  snn.cov.init <- init.snn.cov(fit.cov)
  
  # fit flash for K=1 with the diagonal component, with point exponential prior on L
  fit.cov <- fit.ebcovmf(dat=dat, fl=fit.cov, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1), verbose = 0)$fl
  
  # add one more factor
  K <- fit.cov$n.factors
  fit.cov <- flash.add.greedy(fit.cov, Kmax=1, prior.family = prior.point.laplace(), verbose.lvl = 1)
  
  # add more factors in a greedy manner
  while(fit.cov$n.factors <= Kmax & fit.cov$n.factors==K+1){
    # backfit the added factor
    fit.cov <- flash.backfit(fit.cov, kset=K+1, verbose.lvl = 1)
    fit.cov <- fit.ebcovmf.kset(dat=dat, fl=fit.cov, prior=prior.point.laplace(), kset=K+1)$fl
    
    # initialize the snn fit based on the unconstrained flash fit 
    snn.cov.init <- init.snn.cov(fit.cov, kset=fit.cov$n.factors)
    
    # fit snn with point exponential prior to XXt without considering the diagonal component for now
    fit.cov <- flash.init(dat, var.type = 0) %>% flash.init.factors(EF=snn.cov.init, prior.family=as.prior(ebnm::ebnm_point_exponential, sign = 1)
    ) %>% flash.backfit(verbose.lvl = 0)
    
    # fit snn again with the diagonal component
    fit.cov <- fit.ebcovmf(dat=dat, fl=fit.cov, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1), verbose = 0)$fl
    
    # add one more factor
    K <- fit.cov$n.factors
    fit.cov <- flash.add.greedy(fit.cov, Kmax=1, prior.family = prior.point.laplace(), verbose.lvl = 1) 
  }
  
  # remove the last factor added greedily
  if(fit.cov$n.factors==K+1){
    fit.cov <- flash.remove.factors(fit.cov, kset=fit.cov$n.factors)
  }
  
  return(fit.cov)
} 


### run flash to the covariance matrix by adding factors greedily one by one, before backfitting all existing nonnegative factors, another variant
fit.ebcovmf.snn.v2 <- function(dat, Kmax){
  # fit unconstrained flash for K=1 without considering the diagonal component 
  fit.cov <- flash.init(dat, var.type = 0) %>% flash.add.greedy(Kmax = 1, prior.family = prior.point.laplace(), verbose.lvl = 1
  ) %>% flash.backfit(verbose.lvl = 1)
  
  # initialize the snn fit based on the unconstrained flash fit 
  snn.cov.init <- init.snn.cov(fit.cov)
  
  # fit flash for K=1 with the diagonal component, with point exponential prior on L
  fit.cov <- fit.ebcovmf(dat=dat, fl=fit.cov, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1))$fl
  
  # add one more factor
  K <- fit.cov$n.factors
  fit.cov <- flash.add.greedy(fit.cov, Kmax=1, prior.family = as.prior(ebnm::ebnm_point_exponential, sign = 1), verbose.lvl = 1)
  
  # add one more factor in a greedy manner
  while(fit.cov$n.factors <= Kmax & fit.cov$n.factors==K+1){
    # backfit the added factor
    fit.cov <- flash.backfit(fit.cov, kset=K+1, verbose.lvl = 1)
    fit.cov <- fit.ebcovmf.kset(dat=dat, fl=fit.cov, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1), kset=K+1)$fl
    
    # backfit all factors
    fit.cov <- fit.ebcovmf(dat=dat, fl=fit.cov, prior=as.prior(ebnm::ebnm_point_exponential, sign = 1), verbose = 0)$fl
    
    # add one more factor
    K <- fit.cov$n.factors
    fit.cov <- flash.add.greedy(fit.cov, Kmax=1, prior.family = as.prior(ebnm::ebnm_point_exponential, sign = 1), verbose.lvl = 1) 
  }
  
  # remove the last factor added greedily
  if(fit.cov$n.factors==K+1){
    fit.cov <- flash.remove.factors(fit.cov, kset=fit.cov$n.factors)
  }
  
  return(fit.cov)
} 