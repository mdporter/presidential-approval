## Estimate change points in 538's presidential approval rating

library(tidyverse)
source("R/functions.R")

#-----------------------------------------------------------------------------#
# Polling Data from 538
#-----------------------------------------------------------------------------#
#: load approval data
approval = read_rds("data/approval.rds")

#: Format data
date_start = min(approval$startdate)
data = approval %>% 
  transmute(Y = round(Y), 
            N, 
            L = date2time(startdate, date_start),
            R = date2time(enddate, date_start))


#-----------------------------------------------------------------------------#
# Fit Models to possible change points
#-----------------------------------------------------------------------------#

#: settings
nT = max(data$R)    # observation window [1, nT] (in days)
gap = 7             # minimum days between change points 
lims = c(2, nT-1)   # minimum and maximum possible change points (in days)


#: Fit models and save results
data %>% 
  fit_models(k=0, gap=gap, lims=lims) %>% 
  write_rds("results/k0.rds")

data %>% 
  fit_models(k=1, gap=gap, lims=lims) %>% 
  write_rds("results/k1.rds")

data %>% 
  fit_models(k=2, gap=gap, lims=lims) %>% 
  write_rds("results/k2.rds")

data %>% 
  fit_models(k=3, gap=gap, lims=lims) %>% 
  write_rds("results/k3.rds")

data %>% 
  fit_models(k=4, gap=gap, lims=lims) %>% 
  write_rds("results/k4.rds")

data %>%
 fit_models(k=5, gap=gap, lims=lims) %>%
 write_rds("results/k5.rds")





