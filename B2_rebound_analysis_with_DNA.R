# 06/01/2020
# B2 rebound survival analysis
# with DNA copies and viral_load_2

## package and directory setup
library(tidyverse)
library(survival)
library(survminer)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# 1. read in the previously cleaned and transformed data
dat_log = readRDS("reboundB2_logTrans.rds")

## (somehow log_peak_vl_2 is not in log scale; fix it)
dat_log$log_peak_vl_2 = log(dat_log$log_peak_vl_2, base=10)
## save data again
saveRDS(dat_log, "reboundB2_logTrans.rds")

## also read in DNA load data