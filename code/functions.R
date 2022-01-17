### define the ebnm function for prior (1-pi)*delta_0 + pi*N+(mu, sigma^2), where N+ represents truncated normal on [0, \infty)
ebnm_npmle_truncnorm <- function(x, s, g_init, output=ebnm::output_default()){
  s2 <- s^2
  pi <- g_init$pi[2]
  mu <- g_init$mean[2]
  sigma2 <- g_init$sd[2]^2
  
  ############################## calculate the posterior which is a mixture of delta_0 and truncated normal ##############################
  # calculate mean and variance of posterior corresponding to the truncated normal component in the mixture prior
  wgts <- s2/(s2 + sigma2)
  post_mean <- wgts*mu + (1-wgts)*x
  post_var <- 1/(1/sigma2 + 1/s2)
  post_sd <- pmax(sqrt(post_var), 1e-15)
  
  # calculate the marginal log likelihood for each mixture component
  llik_mat <- matrix(NA, nrow=length(x), ncol=2)
  llik_mat[,1] <- -0.5*(log(s2) + x^2/s2)
  llik_mat[,2] <- -0.5*(log(sigma2+s2) + (x-mu)^2/(sigma2+s2)) + 
    log(pnorm(-post_mean/post_sd, lower.tail = FALSE)) - log(pnorm(-mu/sqrt(sigma2), lower.tail = FALSE))
  llik_norms <- apply(llik_mat, 1, max)
  L_mat <- exp(llik_mat - llik_norms)
  
  # calculate the posterior weight for truncated normal in the posterior
  zeta_mat <- t(t(L_mat)*g_init$pi)
  zeta_mat <- zeta_mat*(1/rowSums(zeta_mat)) 
  zeta <- zeta_mat[,2]
  
  # create the list for return
  ebnm_res <- list()
  
  # return the prior mixture
  if("fitted_g" %in% output)
    ebnm_res$fitted_g <- g_init
  
  # calculate the log marginal likelihood 
  if ("log_likelihood" %in% output) 
    ebnm_res$log_likelihood <- sum(log(L_mat %*% g_init$pi)) + sum(llik_norms) - 0.5*length(x)*log(2*3.141593)
  
  # calculate posterior of theta 
  if ("posterior_mean" %in% output) {
    tmp1 <- truncnorm::etruncnorm(a=0, mean=post_mean, sd=post_sd)
    tmp1[is.nan(tmp1)] <- 0
    tmp1[tmp1 < 0] <- 0
    
    tmp2 <- truncnorm::vtruncnorm(a=0, mean=post_mean, sd=post_sd)
    tmp2[is.nan(tmp2)] <- 0
    tmp2[tmp2 < 0] <- 0
    
    ebnm_res$posterior$mean <- zeta*tmp1
    ebnm_res$posterior$second_moment <- zeta*(tmp1^2 + tmp2)
    ebnm_res$posterior$sd <- sqrt(pmax(0, ebnm_res$posterior$second_moment - ebnm_res$posterior$mean^2))
  }
  
  return(ebnm_res)
}


### define the generalized binary prior (1-pi)*delta_0 + pi*N+(mu, sigma^2), where N+ represents truncated normal on [0, \infty)
### explore if we can improve optimization by spliting it into several subproblems with different range of pi to optimize over
ebnm_binary_general <- function(x, s, g_init, fix_g, output) {
  if (!fix_g) {
    ### objective function to minimize
    opt_fn <- function(par) {
      w <- exp(par[1])/(exp(par[1]) + 1)
      mu <- exp(par[2])
      sigma <- mu/exp(par[3])
      g <- ashr::normalmix(c(1-w, w), c(0, mu), c(0, sigma))
      
      ebnm_res <- ebnm_npmle_truncnorm(x, s, g_init = g, output = "log_likelihood")
      return(-ebnm_res$log_likelihood)
    }
    
    # separately optimize over parameters in different intervals for pi and take the best solution
    pi_list <- c(1e-8, 1e-2, 8e-1)
    wlist <- log(pi_list/(1-pi_list))
    opt_list <- list(NULL)
    val_list <- rep(NA, length(wlist)-1)
    for(k in 1:(length(wlist)-1)){
      opt_list[[k]] <- optim(par=c(wlist[k+1], 0, log(10)), fn = opt_fn, 
                             lower = c(wlist[k], -10, log(10)), upper=c(wlist[k+1], 10, log(100)), method = "L-BFGS-B")
      val_list[k] <- opt_list[[k]]$value
    }
    opt_res <- opt_list[[which.min(val_list)]]
    w <- exp(opt_res$par[1])/(exp(opt_res$par[1]) + 1)
    mu <- exp(opt_res$par[2])
    sigma <- mu/exp(opt_res$par[3])
    g_init <- ashr::normalmix(c(1-w, w), c(0, mu), c(0, sigma))
  }
  
  return(ebnm_npmle_truncnorm(x, s, g_init = g_init, output = output))
}


## initialize the snn fit using flash fit with nonnegative constraints on L
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


### initialize the nonnegative fit to covariance matrix XX' s.t. E[XX'] = LL'+ D based on unconstrained estimate of L
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
fit.ebcovmf <- function(dat, fl, prior=prior.point.laplace(), method="extrapolate", maxiter=500, epsilon=0, verbose=1){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf
  
  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-3) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    kset <- fl$pve > epsilon
    prior.list <- list(NULL)
    prior.list[[1]] <- prior.nonnegative()
    if(sum(kset) > 1){
      for(k in 2:sum(kset)){
        prior.list[[k]] <- prior
      }      
    }
    fl <- flash.init(dat_minuss2) %>% flash.init.factors(EF = lapply(fl$flash.fit$EF, function(x) x[, kset, drop = FALSE]),
          EF2 = lapply(fl$flash.fit$EF2, function(x) x[, kset, drop = FALSE]), prior.family = prior.list) %>% flash.backfit(
          method = method, maxiter = maxiter, verbose.lvl = verbose)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat=dat, fl=fl, s2=s2))
}