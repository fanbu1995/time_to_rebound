# 05/29/2020
# extract viral load data

# 06/03/2020
# per Cliburn's suggestion: 
# calculate AUC of viral load curve pre-ART

# load libraries and set working directory
library(tidyverse)
library(ggplot2)
library(zoo)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# read in the data file (converted to CSV)
VL = read_csv("PVL_B2_02172020_CG.csv")

# create log_viral_load and clean up a bit
VL = VL %>% mutate(log_viral_load = log(viral_load)) %>%
  select(-c(13:14))

saveRDS(VL, "B2_viral_load.rds")

## second-time use: read in the saved data
VL = readRDS("B2_viral_load.rds")

# visualize by animal IDs
ggplot(data=VL, aes(x=week_infection, y=log_viral_load)) + 
  geom_vline(xintercept = c(8,62), size=1) +
  geom_line(aes(color=animal_id)) +
  labs(x="Weeks post infection", 
       y="Log viral load",
       color="Animal ID")+
  scale_x_continuous(breaks = c(0,8,25,50,62,75)) +
  theme_bw(base_size = 14)

# then visualize only the pre-treatment part
ggplot(data=VL%>%filter(week_infection <= 62), 
       aes(x=week_infection, y=log_viral_load)) + 
  geom_vline(xintercept = 8, size=1) +
  geom_line(aes(color=animal_id)) +
  labs(x="Weeks post infection", 
       y="Log viral load",
       color="Animal ID")+
  scale_x_continuous(breaks = c(0,8,25,50)) +
  theme_bw(base_size = 14)

# get "peak viral load" PRE treatment
# by averaging between the largest two values
# and also acquire the log of that averaged value
get_top_mean <- function(x, top=2){
  mean(sort(x,decreasing = T)[1:top], na.rm = T)
}

# here: use "log_10" (to rematch the previous convention...)
peak_VL = VL %>% filter(week_infection <= 62) %>%
  group_by(animal_id) %>%
  summarise(peak_vl_2 = get_top_mean(viral_load)) %>%
  mutate(log_peak_vl_2 = log(peak_vl_2, base=10))

# save it
saveRDS(peak_VL, "peak_viral_load_avg_top2.rds")

# also combine it with the existing log transformed dataset
dat_log = readRDS("reboundB2_logTrans.rds")
dat_log$animal_id = as.character(dat_log$animal_id)
dat_log$group = as.character(dat_log$group)

dat_log$log_peak_vl_2 = peak_VL$peak_vl_2

# and also save this one
saveRDS(dat_log, "reboundB2_logTrans.rds")


# 06/03/2020
# get AUC of viral load

# a function to calculate AUC first
get_AUC <- function(vl, times){
  times = times[-which(is.na(vl))]
  vl = vl[-which(is.na(vl))]
  id = order(times)
  sum(diff(times[id])*zoo::rollmean(vl[id],2))
}

VL_AUC = VL %>% filter(week_infection < 62) %>%
  group_by(animal_id) %>%
  summarise(vl_auc = get_AUC(viral_load, week_infection),
            log_vl_auc = get_AUC(log_viral_load, week_infection))

# save it too
saveRDS(VL_AUC, 'viral_load_AUC.rds')

# also combine this with the previous version of data
dat_log = readRDS("reboundB2_logTrans_withDNA.rds")
dat_log$log_vl_auc = VL_AUC$log_vl_auc

# save it
saveRDS(dat_log, "reboundB2_logTrans_withVLAUC.rds")
