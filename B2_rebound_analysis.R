# 05/28/2020
# first look at the B2 rebound data (cleaned by Yishu Gong)

## package and directory setup
library(tidyverse)
library(survival)
library(survminer)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# read in the data and process/re-organize a little bit
dat = read.csv("reboundp01.csv")
names(dat)

names(dat)[16] = "rebound_time_days_post_ati"

# plot KM curve
simple_KM = survfit(Surv(rebound_time_days_post_ati, observed)~1, data=dat)
plot(simple_KM)
## prettier versions
ggsurvplot(simple_KM, data = dat, 
           censor.shape="|", censor.size = 4,
           size = 1.5, palette = c("#2E9FDF"),
           conf.int = TRUE,
           ggtheme = theme_bw())

ggsurvplot(
  simple_KM, 
  data = dat, 
  size = 1.5,               # change line size
  censor.shape="|",
  censor.size = 4,
  palette = 
    c("#2E9FDF"),           # custom color palettes
  conf.int = TRUE,          # Add confidence interval
  #pval = TRUE,              # Add p-value (to compare two groups...)
  risk.table = TRUE,        # Add risk table
  #risk.table.col = "strata",# Risk table color by groups
  #legend.labs = 
  #  c("Male", "Female"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
