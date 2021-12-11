## Simulation

library(tidyverse)
source("R/functions.R")

#-----------------------------------------------------------------------------#
# Polling Data from 538
#-----------------------------------------------------------------------------#
#-- load approval data
approval = read_rds("data/approval.rds")

#-- Format data
date_start = min(approval$startdate)
data = approval %>% 
  transmute(Y = round(Y), 
            N, 
            L = date2time(startdate, date_start),
            R = date2time(enddate, date_start))

#-----------------------------------------------------------------------------#
# Simulation
#-----------------------------------------------------------------------------#

#: simulation change points and parameters
chgpts = list(
  k0 = NA,
  k1 = 175,
  k2 = c(84, 112), 
  k3 = c(77, 175, 238)
)

alpha = map(chgpts, ~cp_poll(data, .x)$theta$theta)


#: run simulation: generate data and fit all models
run_simulation <- function(data, alpha, chgpts, seed){
  
  #: search space settings
  nT = max(data$R)    # observation window [1, nT] (in days)
  gap = 7             # minimum days between change points 
  lims = c(2, nT-1)   # minimum and maximum possible change points (in days)
  
  #: prior on change points
  p = .4  
  prior_k = tibble(k = 0:5, prior_k = dgeom(k, prob=!!p)) %>% 
    mutate(prior_k = prior_k/sum(prior_k))

  #: Simulate data
  sim = simulate_polls(data, alpha, seed = seed)
  data_sim = sim$data

  #: fit models
  k0 = data_sim %>% 
    fit_models(k=0, gap=gap, lims=lims)
  
  k1 = data_sim %>% 
    fit_models(k=1, gap=gap, lims=lims)
  
  k2 = data_sim %>% 
    fit_models(k=2, gap=gap, lims=lims)
  
  k3 = data_sim %>% 
    fit_models(k=3, gap=gap, lims=lims)
  
  k4 = data_sim %>% 
    fit_models(k=4, gap=gap, lims=lims)
  
  #: results
  B_wide = bind_rows(`0` = k0, `1` = k1, `2` = k2, `3` = k3, `4` = k4, .id="k") %>% 
    mutate(k = as.integer(k)) %>% 
    left_join(prior_k, by="k") %>% 
    group_by(k) %>% mutate(prior_tau = 1/n()) %>% ungroup() %>%  # add prior_tau
    mutate(p = exp(-BIC/2)*prior_k*prior_tau, p = p/sum(p)) %>%  # add posterior prob
    arrange(-p)
  
  
  #: distribution of estimated number of change points
  prob_num_chgpts = B_wide %>% count(k, wt=p)
  
  
  #: performance of true model
  perf_true = B_wide %>% 
    mutate(d1 = abs(t1 - chgpts[1]), 
           d2 = abs(t2 - chgpts[2]), 
           d3 = abs(t3 - chgpts[3])) %>% 
    mutate(rnk = percent_rank(-p)) %>% 
    arrange(d1, d2, d3) %>% 
    slice(1) %>% select(k, p, rnk)
    

  #: save simulation results
  k = sum(!is.na(chgpts))
  save_path = glue::glue("simulation/k{k}/{seed}.rds")
  lst(prob_num_chgpts, perf_true, chgpts, seed) %>% 
    write_rds(save_path)

}


#: run simulation
for(i in 1:100){
  for(k in 0:3){
    print(str_c("starting seed ", i, " and k = ", k))
    seed = 100*i
    k = str_c("k", k)
    run_simulation(data, alpha[[k]], chgpts[[k]], seed = seed)
  }
}




#-----------------------------------------------------------------------------#
# Analyze Simulation Results
#-----------------------------------------------------------------------------#

#: Analyze Results
number_chgpts = tibble()
distr_number_chgpts = list()
proportion_better = tibble()
prob_true = tibble()

for(k in 0:3){
  
  files = list.files(str_c("simulation/k", k), full.names=TRUE)
  results = map(files, read_rds)
  
  # p.correct is the proportion of correctly getting the true number of changepoints
  a = results %>% map_dfr(~.$prob_num_chgpts %>% slice_max(n)) %>%
    summarize(p.correct = mean(k == !!k), n=n()) %>% mutate(k, .before=1)
  number_chgpts = bind_rows(number_chgpts, a)
  
  b = results %>% map_dfr(~.$prob_num_chgpts %>% slice_max(n)) %>%
    count(k) %>% rename(est_k = k) %>% mutate(true_k = !!k, .before=1)
  distr_number_chgpts = bind_rows(distr_number_chgpts, b)
  
  # rnk is the proportion of models that scored better
  c = results %>% map_dfr("perf_true") %>% 
    summarize(pr.better = mean(rnk), n=n()) %>% mutate(k, .before=1)
  proportion_better = bind_rows(proportion_better, c) 
   
  # p is the posterior probability assigned to the true model
  d = results %>% map_dfr("perf_true") %>% 
    summarize(mean.p = mean(p), sd = sd(p), n=n())
  prob_true = bind_rows(prob_true, d)
  
}

number_chgpts
distr_number_chgpts
proportion_better
prob_true






