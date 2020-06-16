# 06/16/2020
# B2 rebound analysis with cell counts and A01 (MHC class)

# as supplement to previous analysis

library(tidyverse)
library(survival)
library(survminer)
library(glmnet)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

dat_log = readRDS("reboundB2_logTrans_CellCounts_A01.rds")

# 1. check out "A01" pos/neg
# 1.1 KM curve stratified by A01 classification
A01_KM = survfit(Surv(rebound_time_days_post_ati, observed)~A01, data=dat_log)
## prettier versions
ggsurvplot(A01_KM, data = dat_log, 
           censor.shape="|", censor.size = 4,
           size = 1.5, palette = c("#E7B800", "#2E9FDF"),
           conf.int = FALSE,
           ggtheme = theme_bw())
