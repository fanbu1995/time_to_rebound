# 05/29/2020
# extract viral load data

# load libraries and set working directory
library(tidyverse)
library(ggplot2)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# read in the data file (converted to CSV)
VL = read_csv("PVL_B2_02172020_CG.csv")


# create log_viral_load
VL = VL %>% mutate(log_viral_load = log(viral_load))

# visualize by animal IDs
ggplot(data=VL, aes(x=week_infection, y=log_viral_load)) + 
  geom_vline(xintercept = c(8,62), size=1) +
  geom_line(aes(color=animal_id)) +
  labs(x="Weeks post infection", 
       y="Log viral load",
       color="Animal ID")+
  scale_x_continuous(breaks = c(0,8,25,50,62,75)) +
  theme_bw(base_size = 14)
  
