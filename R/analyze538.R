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
# Analyze Results
# Note:
#   - it is assumed the same change points are considered for all k. 
#-----------------------------------------------------------------------------#

#: Prior on # change points; geometric distribution (truncated)
p = .4       # geometric parameter
prior_k = tibble(k = 0:4, prior_k = dgeom(k, prob=!!p)) %>% 
  mutate(prior_k = prior_k/sum(prior_k))


#: Load results
tmp = bind_rows(
  `0` = read_rds("results/k0.rds"),
  `1` = read_rds("results/k1.rds"),
  `2` = read_rds("results/k2.rds"),
  `3` = read_rds("results/k3.rds"),
  `4` = read_rds("results/k4.rds"),
  .id = "k"
) %>% 
  mutate(k = as.integer(k)) %>% 
  select(k, t1, t2, t3, t4, BIC)


#: Calculate Probabilities P(model | data)
B_wide = tmp %>%
  left_join(prior_k, by="k") %>%                               # add prior_k
  group_by(k) %>% mutate(prior_tau = 1/n()) %>% ungroup() %>%  # add prior_tau
  mutate(p = exp(-BIC/2)*prior_k*prior_tau, p = p/sum(p)) %>%  # add posterior prob
  arrange(-p)
         

#: number of models
nrow(B_wide)

#: plot-- Pr(k change points)
B_wide %>% 
  count(k, wt=p) %>% 
  ggplot(aes(k, n)) + 
  geom_col() + 
  scale_x_continuous(breaks=0:5, expand=c(0, .1)) + 
  scale_y_continuous(expand=c(.02, 0), breaks=seq(0, 1, by=.1)) + 
  labs(x = "number of change points", y = "probability") + 
  theme_bw()
  

#: Plot-- best fit from each k
data_agg = data %>% 
  rowwise() %>% 
  summarize(t = L:R, len=length(t),  n = N/len, y=Y/len, .groups = "drop") %>% 
  group_by(t) %>% 
  summarize(n = sum(n), y=sum(y), theta=y/n, npolls=n(), .groups="drop") %>% arrange(t) %>% 
  mutate(date = time2date(t, date_start), day=t)


top_models = B_wide %>% mutate(k = as.character(k)) %>% 
  group_by(k) %>% slice_max(p) %>% ungroup() %>% 
  nest(chgpts = matches("t\\d+")) %>% 
  mutate(chgpts = map(chgpts, ~unlist(.x)),
         fit = map(chgpts, ~cp_poll(data, chgpts=.x)),
         theta = map(fit, ~.$theta)) %>% 
  unnest(c(k, theta)) %>% 
  mutate(date = time2date(day, date_start)) 

top_models %>% 
  ggplot(aes(date, theta)) + 
  geom_point(data=data_agg, size=.9) + 
  geom_line(aes(color=k)) + 
  labs(x="date", y="approval") + 
  theme_bw(base_size = 10) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "top") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_colour_brewer(type = "qual", palette="Dark2", name="num of change points") + 
  guides(color = guide_legend(nrow = 1))




#: Get long format to focus on change points
B = B_wide %>% select(-prior_k, -prior_tau) %>% 
  gather(chgpt, time, matches("t\\d+")) %>% 
  filter(!is.na(time)) %>% 
  mutate(date = time2date(time, date_start))


#: Plot-- distribution of change points
B %>% count(date, wt=p, sort=TRUE, name="p") %>% 
  ggplot(aes(date, p)) + geom_area(alpha=.15, color="black") + 
  geom_point() + 
  labs(x="", y="probability") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 


#: Plot-- heatmap of 2 change points
B_wide %>% filter(!is.na(t1), !is.na(t2)) %>% 
  count(t1, t2, wt=p, name="p") %>% 
  mutate(p = p/sum(p)) %>% 
  mutate_at(vars(-p), ~time2date(.x, date_start)) %>% 
  ggplot(aes(t1, t2, fill=p)) + 
  geom_tile() + 
  scale_fill_viridis_c(option="D", direction=1) + 
  geom_point(data=. %>% top_n(n=1, wt=p), color="red", size=3) +
  theme_bw() + labs(x="1st change point", y="2nd change point") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + 
  scale_y_date(date_labels = "%b", date_breaks = "1 month")



# Other plots not used in paper

#: Colored plot of change points
B %>% 
  group_by(chgpt, date) %>% summarize(p = sum(p), .groups = "drop")  %>% 
  mutate(chgpt = recode(chgpt, t1 = "1st", t2 = "2nd", t3 = "3rd", t4 = "4th")) %>% 
  ggplot(aes(date, p)) + 
  geom_area(aes(fill=chgpt, color=chgpt), alpha=.15, 
            position="identity") + 
  geom_line(data=. %>% group_by(date) %>% summarize(p = sum(p)),
            size=1.05) + 
  labs(x="", y="probability") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 

#: Pr(t_k = t)
B %>% 
  count(chgpt, date, wt=p, name="p") %>% 
  ggplot(aes(date, p)) + geom_area(alpha=.15, color="black") + geom_point() + 
  facet_wrap(~chgpt, ncol=1) + 
  labs(x="", y="probability", title="Pr(t_k = t)") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 


#: Pr(t_k = t, k)
B %>% 
  count(k, date, wt = p, name="p") %>% 
  ggplot(aes(date, p)) + geom_area(alpha=.15, color="black") + 
  geom_point() + 
  facet_wrap(~k, ncol=1) + 
  labs(x="", y="probability", title="Pr(t_k = t, k)") + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw(base_size = 10) + theme(axis.title.x = element_blank()) 
  

  







