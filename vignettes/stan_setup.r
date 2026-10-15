## Setup for stan models via brms
## Defining negative binomial with multiple trials
suppressPackageStartupMessages({
    require(brms)
})

neg_binomial_2_trials <- custom_family(
    "neg_binomial_2_trials", 
    dpars = c("mu", "phi"),
    links = c("log", "log"), 
    lb = c(0, 0),
    type = "int", 
    vars = "trials[n]"
)

neg_binomial_2_trials_phi <- custom_family(
  "neg_binomial_2_trials_phi", 
  dpars = c("mu", "phi"),
  links = c("log", "log"), 
  lb = c(0, 0),
  type = "int", 
  vars = "trials[n]"
)

neg_binomial_2_phi_nu <- custom_family(
  "neg_binomial_2_phi_nu", 
  dpars = c("mu", "phi", "nu"),
  links = c("log", "log", "identity"), 
  lb = c(0, 0, 0),
  type = "int", 
  vars = "trials[n]"
)

neg_binomial_2_nu <- custom_family(
    "neg_binomial_2_nu", 
    dpars = c("mu", "phi", "nu"),
    links = c("log", "log", "identity"), 
    lb = c(0, 0, 0),
    type = "int", 
    vars = "trials[n]"
)

stan_funs <- "
    real neg_binomial_2_trials_lpmf(int y, real mu, real phi, int T) {
        return neg_binomial_2_lpmf(y | mu * T, phi * T);
    }
     real neg_binomial_2_trials_lccdf(int y, real mu, real phi, int T) {
        return neg_binomial_2_lccdf(y | mu * T, phi * T);
     }
     real neg_binomial_2_trials_lcdf(int y, real mu, real phi, int T) {
        return neg_binomial_2_lcdf(y | mu * T, phi * T);
    }
    int neg_binomial_2_trials_rng(real mu, real phi, int T) {
        return neg_binomial_2_rng(mu * T, phi * T);
    }
     real neg_binomial_2_trials_phi_lpmf(int y, real mu, real phi, int T) {
        return neg_binomial_2_lpmf(y | mu * T, phi * T);
    }
    int neg_binomial_2_trials_phi_rng(real mu, real phi, int T) {
        return neg_binomial_2_rng(mu * T, phi * T);
    }
    real neg_binomial_2_nu_lpmf(int y, real mu, real phi, real nu, int T) {
        return neg_binomial_2_lpmf(y | mu * T, phi * T * (mu^(2*nu)));
    }
    int neg_binomial_2_nu_rng(real mu, real phi, real nu, int T) {
        return neg_binomial_2_rng(mu * T, phi * T * (mu^(2*nu)));
    }
    real neg_binomial_2_phi_nu_lccdf(int y, real mu, real phi, real nu, int T) {
        return neg_binomial_2_lccdf(y | mu * T, phi * T * (mu^(2*nu)));
    }
    real neg_binomial_2_phi_nu_lcdf(int y, real mu, real phi, real nu, int T) {
        return neg_binomial_2_lcdf(y | mu * T, phi * T * (mu^(2*nu)));
    }
    real neg_binomial_2_phi_nu_lpmf(int y, real mu, real phi, real nu, int T) {
        return neg_binomial_2_lpmf(y | mu * T, phi * T * (mu^(2*nu)));
    }
    int neg_binomial_2_phi_nu_rng(real mu, real phi, real nu, int T) {
        return neg_binomial_2_rng(mu * T, phi * T * (mu^(2*nu)));
    }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")

## For prediction:
#expose_functions(model, vectorize = TRUE)

log_lik_neg_binomial_2_trials <- function(i, draws){
    mu <- draws$dpars$mu[, i]
    phi <- draws$dpars$phi
    trials <- draws$data$trials[i]
    y <- draws$data$Y[i]
    neg_binomial_2_trials_lpmf(y, mu, phi, trials)
}

posterior_predict_neg_binomial_2_trials <- function(i, draws, ...) {
  mu <- pmax(draws$dpars$mu[, i], 1E-10)
  if(any(mu>100)) cat("i: ",i,"\n")
  mu <- pmin(mu, 100)
  phi <- pmax(draws$dpars$phi, 1E-10)
  trials <- pmax(draws$data$trials[i], 1)
  
  return(neg_binomial_2_trials_rng(mu, phi, trials))
}

posterior_epred_neg_binomial_2_trials <- function(prep) {
  return(brms:::posterior_epred_negbinomial2(prep))
}

log_lik_neg_binomial_2_trials_phi <- function(i, draws){
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi[, i]
  trials <- draws$data$trials[i]
  y <- draws$data$Y[i]
  neg_binomial_2_trials_phi_lpmf(y, mu, phi, trials)
}

posterior_predict_neg_binomial_2_trials_phi <- function(i, draws, ...) {
  mu <- pmax(draws$dpars$mu[, i], 1E-10)
  phi <- pmax(draws$dpars$phi[, i], 1E-10)
  trials <- pmax(draws$data$trials[i], 1)
  return(neg_binomial_2_trials_phi_rng(mu, phi, trials))
}

posterior_epred_neg_binomial_2_trials_phi <- function(prep) {
  return(brms:::posterior_epred_negbinomial2(prep))
}

log_lik_neg_binomial_2_nu <- function(i, draws){
    mu <- draws$dpars$mu[, i]
    phi <- draws$dpars$phi
    nu <- draws$dpars$nu    
    trials <- draws$data$trials[i]
    y <- draws$data$Y[i]
    neg_binomial_2_nu_lpmf(y, mu, phi, nu, trials)
}

posterior_predict_neg_binomial_2_nu <- function(i, draws, ...) {
  mu <- pmax(draws$dpars$mu[, i], 1E-10)
  mu <- pmax(draws$dpars$mu[, i], 1E-10)
  if(any(mu>100)) cat("i: ",i,"\n")
  mu <- pmin(mu, 100)
  phi <- pmax(draws$dpars$phi, 1E-10)
  nu <- pmax(draws$dpars$nu, 1E-10)  
  trials <- pmax(draws$data$trials[i], 1)
  return(neg_binomial_2_nu_rng(mu, phi, nu, trials))
}

posterior_epred_neg_binomial_2_nu <- function(prep) {
  return(brms:::posterior_epred_negbinomial2(prep))
}

log_lik_neg_binomial_2_phi_nu <- function(i, draws){
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi[, i]
  nu <- draws$dpars$nu    
  trials <- draws$data$trials[i]
  y <- draws$data$Y[i]
  neg_binomial_2_phi_nu_lpmf(y, mu, phi, nu, trials)
}

posterior_predict_neg_binomial_2_phi_nu <- function(i, draws, ...) {
  
  mu <- pmax(draws$dpars$mu[, i], 1E-10)
  if(any(mu>100)) cat("i: ",i,"\n")
  mu <- pmin(mu, 100)
  phi <- pmax(draws$dpars$phi[, i], 1E-10)
  nu <- pmax(draws$dpars$nu, 1E-10)  
  trials <- pmax(draws$data$trials[i], 1)
  return(neg_binomial_2_phi_nu_rng(mu, phi, nu, trials))
}

posterior_epred_neg_binomial_2_phi_nu <- function(prep) {
  return(brms:::posterior_epred_negbinomial2(prep))
}
#loglik <- log_lik(model)
#looic <- loo(model)
#save(looic, loglik, file='predictions.rdata')
