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


### apply unconstrained flash to covariance matrix XX' s.t. E[XX'] = LL'+ D, where D = sigma2*I
fit.ebcovmf <- function(dat, fl, prior=prior.point.laplace()){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf
  
  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    fl <- flash.init(dat_minuss2) %>% flash.init.factors(EF = fl$flash.fit$EF, EF2 = fl$flash.fit$EF2, prior.family = prior) %>% flash.backfit(verbose.lvl = 1)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat=dat, fl=fl, s2=s2))
}


## initialize the nonnegative fit to covariance matrix XX' s.t. E[XX'] = LL'+ D based on unconstrained estimate of L
init.snn.cov <- function(fl) {
  LL <- fl$flash.fit$EF[[1]]
  FF <- fl$flash.fit$EF[[2]]
  
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  
  to.keep <- (colSums(LL) > .Machine$double.eps) & (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  
  return(list(LL, FF))
}