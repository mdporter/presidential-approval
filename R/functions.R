## Functions for estimating change points in polls

# Utility Functions
#=============================================================================#
# Convert date to integer time and time back to date
#=============================================================================#
date2time <- function(date, start_date, ...){
  start_date = as.Date(start_date)
  as.integer(as.Date(date, ...) - start_date + 1)
}
time2date <- function(time, start_date){
  start_date = as.Date(start_date)
  time + start_date - 1
}



# joinpoint_matrix()
#=============================================================================#
# Model matrix for join-point regression
#
# The design matrix is X[i,j] = (x_i - t_j)*(x > t_j). Where t_j is (j-1)th change 
#   point. The first column is for slope (essentially change at time 1). 
# 
# Inputs:
#   x: vector of times. 
#   chgpts: values of potential change points. Corresponds to the start of the 
#           change. So a change at t=2 means the deviation will not be seen until 
#           time 3. If chgpts is NA, then no change points are used.
#   intercept: should column of 1's be added
#
# Outputs:
#   Model/Design matrix
#
# Notes:
#   - model is b0 + b1*x + sum_j b_j (x-t_j)*(x>t_j)
#   - no change point at max(x) since there is no data to estimate
#=============================================================================#
joinpoint_matrix <- function(x, chgpts = seq(min(x)+1, max(x)-1, by=1),
                             intercept = TRUE) {
  
  p = length(chgpts)
  if(p == 1 && is.na(chgpts)) p = 0       # if chgpts=NA, treat as no change point
  X = matrix(0, length(x), p+1)
  X[,1] = x                               # first columns is for slope
  if(p > 0){
    for(j in 1:p) {                       # other columns are for change points
      tau = chgpts[j]
      X[, j+1]  = (x - tau)*(x>tau)
    }
    colnames(X) = c("slope", paste0("chgpt=", chgpts))
  } else colnames(X) = "slope" 
  if(intercept) X = cbind("intercept"=1, X)
  attr(X, "chgpts") = chgpts
  return(X)
}


# cp_poll()
#=============================================================================#
# Linear joinpoint model (i.e., trend filtering) for interval censored Poisson data
#
# Fits a (log) linear joinpoint model at specific change points
# 
# Inputs:
#   data: data frame with columns named Y, N, L, R
#     Y is number of successes
#     N is poll size
#     L,R are start and end days of poll (must be integers)
#   chgpts: values of potential change points. Corresponds to the start of the 
#           change. So a change at time 2 means the deviation will not be seen until 
#           time 3. If chgpts = NA, then linear fit with no change points
#   maxit: maximum number of iterations
#   tol: convergence tolerance for the mean log-likelihood ratio gain. 
#
# Outputs:
#   list with objects
#     logL: log likelihood at convergence
#     BIC: -2logL + p*log(n)
#     pars: estimated parameters
#     niters: number of iterations
#     logL0: log likelihood at initialization
#     theta: tibble of day and estimated rate
#
# Notes:
#   - initializes theta as a constant
#   - requires dplyr (> 1.0.1)
#
# Todo/consider:
#   - let theta.start be input (warm starts)
#=============================================================================#
cp_poll <- function(data, chgpts = c(100, 200), 
                    maxit = 25, tol=1E-7) {
  
  #-- Checks
  if( !all(sort(colnames(data)) == sort(c("Y","N","L","R"))) ) 
    stop('data must have column names: Y, N, L, R')
  
  chgpts = chgpts[!is.na(chgpts)] # remove all NA's
  
  if(length(chgpts) > 0) {
    if(min(chgpts) < (min(data$L)+1)) stop('chgpts must be after than ealiest poll')
    if(max(chgpts) > (max(data$R)-1)) stop('chgpts must be before the latest poll')
  }
  
  #-- Set-up
  nT = max(data$R)    # maximum end date
  npolls = nrow(data) # number of polls
  data_EM = data %>% 
    mutate(i = 1:npolls) %>% # add poll index
    group_by(i) %>%          # group by poll index
    summarize(t = L:R, len=length(t),  n = N/len, y=Y, .groups = "keep")
  
  #-- Model Matrices  
  X = joinpoint_matrix(data_EM$t, chgpts = chgpts, intercept = TRUE)
  X_full = joinpoint_matrix(1:nT, chgpts = chgpts, intercept = TRUE)
  
  #-- Initialize theta (constant)
  theta = with(data, sum(Y) / sum(N))
  data_EM$theta = theta
  logL0 = with(data, sum(dpois(Y, lambda=N*theta, log=TRUE)))
  
  #-- EM  
  logL = numeric(maxit)
  Pars = matrix(0, maxit, ncol(X), dimnames = dimnames(X))
  converged = FALSE
  for(i in 1:maxit) {
    
    #-- E-step
    data_EM = data_EM %>% 
      mutate(p = theta*n, p = p/sum(p), # grouped mutate
             z = y*p, 
             z = round(z))
    
    #-- M-step
    fit = glm.fit(x=X, y=data_EM$z, 
                  offset = log(data_EM$n),
                  family=poisson())
    theta = exp(X_full %*% matrix(fit$coefficients)) %>% as.numeric()
    data_EM$theta = theta[data_EM$t]
    
    #-- Keep
    Pars[i,] = fit$coefficients
    logL[i] = data_EM %>%
      summarize(mu = sum(theta*n), y=y[1], .groups="drop") %>% # grouped summarize
      summarize(logL = sum( dpois(y, mu, log=TRUE))) %>% pull()
    
    #-- Convergence check
    if(i>1 && logL[i] - logL[i-1] <  npolls*tol) {
      converged = TRUE
      break
    }
  }
  
  #-- Output
  if(!converged) warning("did not converge; increase maxit")
  
  list(logL = logL[i], 
       BIC = -2*logL[i] + ncol(X)*log(npolls),
       pars = Pars[i,], 
       niters = i, 
       logL0 = logL0,
       theta = tibble(day=seq_along(theta), theta)
  )
}


# search_space()
#=============================================================================#
# Make set of candidate change points
#
# Inputs:
#   k: number of change points {0, 1, ...}
#   gap: minimum time between change points
#   lims: minimum and maximum values of the change points
#
# Outputs:
#   tibble of k change points
# Notes:
#   requires tibble, dplyr, purrr
#=============================================================================#
search_space <- function(k=3, gap=7, lims=c(2, 99)){
  if(any(k < 0)) stop("k must be an integer >= 0")
  x = seq(max(gap, lims[1]), lims[2], by=gap)
  get_combos <- function(m, x) {
    if(m == 0) return(tibble(t1=NA))
    combn(x, m) %>% t() %>% 
      as_tibble(.name_repair = "minimal") %>% 
      setNames(str_c("t", 1:length(.)))
  }
  map_dfr(k, ~get_combos(., x))
}

# fit_models()
#=============================================================================#
# Creates search space, fits models, and records BIC score.
# runs search_space() and cp_poll() 
#
# Inputs:
#   see arguments for search_space() and cp_poll()
#
# Outputs:
#   tibble of change points and their BIC scores
#
# Notes:
#   requires tibble, dplyr, purrr, progress (for progress bar)
#=============================================================================#
fit_models <- function(data, k, gap, lims){
  tau = search_space(k=k, gap=gap, lims=lims) 
  BIC = numeric(nrow(tau))
  pb = progress::progress_bar$new(total = nrow(tau), 
                                  format="[:bar] :elapsed time; :eta remaining")
  for(i in 1:nrow(tau)){
    pb$tick()
    BIC[i] = cp_poll(data, tau[i,])$BIC
  }
  mutate(tau, BIC)
}

# ## requires: progressr, and (optional) furrr
# fit_models_parallel <- function(data, k, gap, lims, parallel = FALSE){
#   tau = search_space(k=k, gap=gap, lims=lims) %>% split(1:nrow(.))
#   progressr::with_progress({
#     p <- progressor(steps = length(tau))
#     if(parallel) BIC <- furrr::future_map_dbl(tau, ~{p(); cp_poll(data, unlist(.x))$BIC})
#     else BIC <- map_dbl(tau, ~{p(); cp_poll(data, unlist(.x))$BIC})
#   })
#   bind_rows(tau) %>% mutate(BIC)
# }
# 
# library(progressr)
# library(furrr)
# 
# availableCores()
# plan(multisession, workers = 8)
# handlers(
#   handler_progress(
#     format   = "Number of models: :total | (:message) [:bar] :percent | Total Time: :elapsed | Remaining: :eta",
#     width    = 80,
#     complete = "+",
#     clear    = FALSE
#   )
# )


# simulate_polls()
#=============================================================================#
# Simulate aggregated polling data
#
# Inputs:
#   data: original observed data
#   chgpts: vector of change points (days)
#   beta: coefficients for support rate 
#    vector of at least length 2. First two elements are intercept and slope.
#    any other elements are slope change at change points
#   nT: maximum number of days under observation
#   seed: random seed. If NULL, then doesn't set seed.
#
# Outputs:
#   a list with elements
#   - data: simulated data (data frame)
#   - alpha: the alpha values per day (data frame)
#   - seed: random seed used to generate data
#
# Details:
#   Uses the same number and date ranges of original poll data. 
#   The number of supporters is Poisson with mean sum alpha*n and number of 
#     non-supporters is Poisson with mean of sum (1-alpha)*n,
#     where alpha = exp(X %*% beta) is the support probability and n is the 
#     number of respondents (assumed to be equally distributed across all 
#     survey days)
#=============================================================================#
simulate_polls <- function(data, chgpts, beta, nT = max(data$R), seed=NULL){

  #-- Checks
  if(length(chgpts) != length(beta) - 2){
    stop("beta must have length(chgpts) + 2 elements")
  }
    
  #-- model matrix
  X_full = joinpoint_matrix(1:nT, chgpts = chgpts, intercept = TRUE)
  
  #-- Make alpha (support probability)
  beta = as.matrix(beta)
  eta = as.numeric(X_full %*% beta)
  alpha = exp(eta)
  if(any(alpha > 1)) stop("alpha shouldn't be greater than one. Change beta.")
  
  #-- Simulate data
  if(!is.null(seed)) set.seed(seed)
  
  data_sim_agg = data %>%
    #: generate new N (sample size)
    mutate(N = rpois(n(), lambda=N)) %>% 
    #: group by poll 
    mutate(id = row_number()) %>% 
    group_by(id, N, L, R) %>% 
    #: add alpha and calculate lambda
    mutate(day = list(L:R)) %>% unnest(day) %>% 
    mutate(alpha = alpha[day], n = N/n()) %>% 
    summarize(lambda = sum(n*alpha), .groups = "drop") %>% 
    #: generate new Y (success counts)
    mutate(Y = rbinom(n(), size=N, prob=lambda/N)) %>% 
    #: select rows
    select(Y, N, L, R)

  list(data = data_sim_agg, 
       alpha = tibble(day = seq_along(alpha), alpha), 
       seed = seed)
}
  
  

