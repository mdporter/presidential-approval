#=============================================================================#
# Retrieve and process polling data from 538
# https://projects.fivethirtyeight.com/trump-approval-ratings/
# source: https://projects.fivethirtyeight.com/trump-approval-data/approval_polllist.csv
#
# Notes: 
#   Using A graded polls only (A+, A, A-)
#   date range: 2020/01/01 - 2020/12/03
#   subgroup: All polls
#   population: using all (adults, likely voters, registered voters)
# Outputs:
#   data/approval.rds 
#=============================================================================#
library(readr)
library(dplyr)


#-- Load 538 approval polls
polls = read_csv("https://projects.fivethirtyeight.com/trump-approval-data/approval_polllist.csv") %>% 
  mutate_at(vars(startdate, enddate), as.Date, format = "%m/%d/%Y")

#-- Settings
poll_grade = c("A+","A","A-")                     # use A graded polls only
date_rng = as.Date(c("2020-01-01", "2020-12-03")) # 2020, up to election day
sub_grp = "All polls"                             # {All polls, Adults, Voters}

#-- approval ratings
approval = polls %>% 
  # Select relevant polls
  filter(president == "Donald Trump",  
         subgroup == sub_grp,          
         grade %in% poll_grade,        
         startdate >= date_rng[1], enddate <= date_rng[2] 
         ) %>% 
  # select columns
  transmute(startdate, enddate, 
            p = adjusted_approve / 100,
            Y = p * samplesize,
            N = samplesize,
            grade,
            poll_id) %>% 
  distinct() # ensure no duplicates


#-- save data
write_rds(approval, "data/approval.rds")
