## Estimate change points in 538 presidential approval rating

library(tidyverse)
source("R/functions.R")

#-----------------------------------------------------------------------------#
# Polling Data from 538
# https://projects.fivethirtyeight.com/trump-approval-ratings/
#-----------------------------------------------------------------------------#
#-- load approval data
approval = load_polls(poll_grade = c("A+","A","A-"), # use A graded polls only)
                      date_rng = c(as.Date("2020-01-01"), as.Date("2020-12-07")))

# write_rds(approval, "data/approval.rds")

#-- Format data
date_start = min(approval$startdate)
data = approval %>% 
  transmute(Y = round(Y), # round so dpois() doesn't give warnings
            N, 
            L = date2time(startdate, date_start),
            R = date2time(enddate, date_start))


#-----------------------------------------------------------------------------#
# Search 0-5 change points
#-----------------------------------------------------------------------------#
nT = max(data$R) # observation window [1, nT] (in days)
gap = 7
lims = c(2, nT-1)

search_space(k=0, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k0.rds")

search_space(k=1, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k1.rds")

search_space(k=2, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k2.rds")

search_space(k=3, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k3.rds")

search_space(k=4, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k4.rds")

search_space(k=5, gap=gap, lims=lims) %>% 
  mutate(BIC = pmap_dbl(., ~cp_poll(data, .x)$BIC)) %>% 
  write_rds("results/k5.rds")

# #-----------------------------------------------------------------------------#
# # Search 0-5 Change points using furrr (and futures) package
# #-----------------------------------------------------------------------------#
# library(furrr)
# library(progressr)
#
# #-- Set up 
# availableCores()
# plan(multisession, workers = 10)
# handlers(
#   handler_progress(
#     format   = "Number of models: :total | (:message) [:bar] :percent | Total Time: :elapsed | Remaining: :eta",
#     width    = 60,
#     complete = "+",
#     clear    = FALSE
#   )
# )
# 
# #-- Search Space
# tau = search_space(k=5, gap=gap, lims=lims)
# 
# #-- Run 
# with_progress({
#   p <- progressor(steps = nrow(tau))
#   BIC <- future_pmap_dbl(tau, ~{p(); cp_poll(data, .x)$BIC})
# })
# A = mutate(tau, BIC)



#-----------------------------------------------------------------------------#
# Analyze Results
# Note:
#   - it is assumed the same change points are considered for all k. 
# TODO: Estimate BIC values for all possible change points
#-----------------------------------------------------------------------------#

#-- Prior on # change points; geometric distribution
p = .4  # geometric parameter
prior_k = tibble(k = 0:5, prior_k = dgeom(k, prob=!!p)) %>% 
  mutate(prior_k = prior_k/sum(prior_k))

g = ggplot(prior_k, aes(k, prior_k)) + scale_x_continuous(breaks=0:5)
g + geom_col()  + labs(y = "Pr(K=k)")  # pmf
g + stat_ecdf() + labs(y = "Pr(K<=k)") # cdf


#-- Load Results
get_results <- function(path="k3.rds", dir="results") {
  d = read_rds(file.path(dir, path)) 
  #- find number of change points
  k = d %>% select(matches("t\\d+")) %>% mutate_all(~!is.na(.x)) %>% rowSums()
  #- return long format
  d %>% mutate(k, model = str_c(k, row_number(), sep="_")) %>% 
      gather(chgpt, time, -BIC, -k, -model)
}

tmp = bind_rows(
  get_results("k0.rds"),
  get_results("k1.rds"),
  get_results("k2.rds"),
  get_results("k3.rds"),
  get_results("k4.rds"),
  get_results("k5.rds")
)

#-- Calculate Probabilities P(model | data)
B_wide = tmp %>% spread(chgpt, time) %>% 
  left_join(prior_k, by="k") %>%                               # add prior_k
  group_by(k) %>% mutate(prior_tau = 1/n()) %>% ungroup() %>%  # add prior_tau
  mutate(p = exp(-BIC/2)*prior_k*prior_tau, p = p/sum(p)) %>%  # add posterior prob
  arrange(-p)
         

#-- number of models
nrow(B_wide)
n_distinct(A$model)

#-- Pr(k change points)
B_wide %>% 
  count(k, wt=p) %>% 
  ggplot(aes(k, n)) + 
  geom_col() + 
  scale_x_continuous(breaks=0:5, expand=c(0, .1)) + 
  scale_y_continuous(expand=c(.02, 0), breaks=seq(0, 1, by=.1)) + 
  labs(x = "number of change points", y = "probability") + 
  theme_bw()
  


#-- Get long format to focus on change points
B = B_wide %>% select(-prior_k, -prior_tau) %>% 
  gather(chgpt, time, matches("t\\d+")) %>% 
  filter(!is.na(time)) %>% 
  mutate(date = time2date(time, date_start))


#-- Marginal Probabilities of change points
B %>% 
  count(chgpt, date, wt=p, name="p") %>% 
  spread(chgpt, p, fill=0) %>% 
  arrange(-t1)

B %>% 
  count(date, wt=p, sort=TRUE, name="p") %>% 
  arrange(-p)


#-- Plot: probability of change point
B %>% count(date, wt=p, sort=TRUE, name="p") %>% 
  ggplot(aes(date, p)) + geom_area(alpha=.15, color="black") + 
  geom_point() + 
  labs(x="", y="probability", title="Probability of change point") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 

B %>% count(chgpt, date, wt=p, name="p") %>% 
  ggplot(aes(date, p)) + 
  geom_area(aes(fill=chgpt, color=chgpt), alpha=.15, 
            position="identity") + 
  geom_line(data=. %>% group_by(date) %>% summarize(p = sum(p)),
            size=1.05) + 
  labs(x="", y="probability", title="Probability of nth change point") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 

# Pr(t_k = t)
B %>% 
  count(chgpt, date, wt=p, name="p") %>% 
  # group_by(chgpt) %>% mutate(p = p/sum(p)) %>% # Pr(t_k = t | K >= k)
  ggplot(aes(date, p)) + geom_area(alpha=.15, color="black") + geom_point() + 
  facet_wrap(~chgpt, ncol=1) + 
  labs(x="", y="probability", title="Probability of change point") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 

  
#-- Plot heatmap of 2-way change points

# conditional on k=2
B_wide %>% filter(k == 2) %>% 
  count(t1, t2, wt=p, name="p") %>% 
  mutate(p = p/sum(p)) %>% 
  mutate_at(vars(-p), ~time2date(.x, date_start)) %>% 
  ggplot(aes(t1, t2, fill=p)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  geom_point(data=. %>% top_n(n=1, wt=p), color="red", size=3) +
  theme_bw() + labs(x="1st change point", y="2nd change point") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + 
  scale_y_date(date_labels = "%b", date_breaks = "1 month")

# unconditional
B_wide %>% filter(!is.na(t1), !is.na(t2)) %>% 
  count(t1, t2, wt=p, name="p") %>% 
  mutate(p = p/sum(p)) %>% 
  mutate_at(vars(-p), ~time2date(.x, date_start)) %>% 
  ggplot(aes(t1, t2, fill=p)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  geom_point(data=. %>% top_n(n=1, wt=p), color="red", size=3) +
  theme_bw() + labs(x="1st change point", y="2nd change point") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + 
  scale_y_date(date_labels = "%b", date_breaks = "1 month")


#-- Plot: best fit from each k

#- aggregated (naively) data
data_agg = data %>% 
  rowwise() %>% 
  summarize(t = L:R, len=length(t),  n = N/len, y=Y/len, .groups = "drop") %>% 
  group_by(t) %>% 
  summarize(n = sum(n), y=sum(y), theta=y/n, npolls=n(), .groups="drop") %>% arrange(t) %>% 
  mutate(date = time2date(t, date_start))


B_wide %>% mutate(k = as.character(k)) %>% 
  group_by(k) %>% slice_max(p) %>% ungroup() %>% 
  nest(chgpts = matches("t\\d+")) %>% 
  mutate(chgpts = map(chgpts, ~unlist(.x)),
         fit = map(chgpts, ~cp_poll(data, chgpts=.x)),
         theta = map(fit, ~.$theta)) %>% 
  unnest(c(k, theta)) %>% 
  mutate(date = time2date(day, date_start)) %>% 
  ggplot(aes(date, theta)) + 
  geom_point(data=data_agg, size=.9) + 
  geom_line(aes(color=k)) + 
  labs(x="date", y="probability") + 
  theme_bw(base_size = 10) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "top") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_colour_brewer(type = "qual", palette="Dark2", name="# chg pts") + 
  guides(color = guide_legend(nrow = 1))
  
  
  




