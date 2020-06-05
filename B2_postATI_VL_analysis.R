# 06/05/2020
# look a bit more at post ATI viral levels
# and see what factors influence them


# load required packages and set working directory
library(tidyverse)
library(ggplot2)
#library(zoo)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")


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

# visualize only the part after ATI
VL_post = VL %>% filter(days_post_ati >= 0)

# load rebound time data as well (mark up rebound time in plots)
dat_log = readRDS("reboundB2_logTrans_withVLAUC.rds")
rebound_times = dat_log %>% 
  select(animal_id, rebound_time_days_post_ati)
rebound_times$rebound_time_days_post_ati[rebound_times$animal_id=="RQc19"] = NA

ggplot(data=VL_post, aes(x=days_post_ati, y=log_viral_load)) +
  geom_hline(yintercept = log(60, base=10), 
             size = 0.5, linetype = "dashed") +
  geom_vline(aes(xintercept = rebound_time_days_post_ati),
             data = rebound_times, 
             size = 0.5, linetype = "dashed") +
  geom_line(aes(color=animal_id)) +
  geom_point(aes(color=animal_id), size=1) +
  labs(x="Days post ATI", 
       y="Viral load (log10)",
       color="Animal ID")+
  #scale_x_continuous(breaks = c(0,8,25,50,62,75)) +
  scale_y_continuous(limits = c(0,8)) +
  theme_bw(base_size = 14)+
  facet_wrap(~animal_id, ncol = 5)
