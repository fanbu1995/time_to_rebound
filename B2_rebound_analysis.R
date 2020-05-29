# 05/28/2020
# first look at the B2 rebound data (cleaned by Yishu Gong)

## package and directory setup
library(tidyverse)
library(survival)
library(survminer)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# 1. read in the data and process/re-organize a little bit
dat = read.csv("reboundp01.csv")
names(dat)

## (clean up the name of rebound_time)
names(dat)[16] = "rebound_time_days_post_ati"

## (remove repeated column for "point_ic50_8_weekspost_ATI")
dat = dat %>% select(-point_ic50_8_weekspost_ATI.1)

## (move the columns on ART interruption, rebound time, and observed indicator
## to column 3:5, so that columns 6:30 are the covariates)
dat = dat %>% select(1:2, 15:17, 3:14, 18:30)

## (change the gp41 and gp120 columns naming conventions for simplicity)
names(dat)[c(10:17)] = c("peak_gp41", "log_peak_gp41",
                         "gp41_treat", "log_gp41_treat",
                         "peak_gp120", "log_peak_gp120",
                         "gp120_treat", "log_gp120_treat")

## (change the "SHIV.RNAgag_copies_permillion_CD4_RNA_copies" column names for simplicity)
names(dat)[18:22] = c('RNA_copies_blood_8', 'RNA_copies_blood_56',
                      'RNA_copies_LN_8', 'RNA_copies_LN_56',
                      'RNA_copies_RB_56') # RB=rectal_biopsy



# 1.5 examine covariate distributions and think about transformations
## Histograms
for(n in names(dat)[6:30]){
  hist(dat[,n], main=n)
}

### Findings:
# - For viral loads, antibodies, the log scale makes much more sense 
#   (otherwise range too large)
# - RNA copies per million: range large (except for LN 56 weeks), 
#   maybe need to log transform too
# - point_ic50: also large range (except for 8 weeks post ATI),
#   maybe need to log transform too
# - pos_auc doesn't need log transformation?? (seems to be between 0 and 1)

## create log scale measurements for:
# - SHIV RNA copies per million CD4 RNA
# - point IC 150
dat = dat %>% 
  mutate_at(18:26, list(log=log)) %>% 
  rename_at(vars( contains( "_log") ), list( ~paste("log", gsub("_log", "", .), sep = "_") ))

## also: create a separate dataset with all the log transformations done
dat_log = dat %>% select(c(1:5,7,9,11,13,15,17,31:39,27:30))

## (save it for usage)
saveRDS(dat_log,"reboundB2_logTrans.rds")


# 2. plot KM curve
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

# 3. examine correlation/relationship between covariates

## (Below: some customized functions for "pairs")
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## 3.1 viral loads and antibody concentration, in log scale
# pairs(dat[,c(7,9,11,13,15,17)], 
#       lower.panel = panel.cor, diag.panel = panel.hist)
pairs(dat_log %>% select(log_peak_vl:log_gp120_treat), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
# - GP41 concentration at peak and at treatment PERFECTLY correlated
# - peak viral load highly correlated with:
#   - GP120 at treatment
#   - GP120 at peak
# - GP120 concentration at peak and at treatment also highly correlated

## 3.2 RNA copies
pairs(dat_log %>% select(log_RNA_copies_blood_8:log_RNA_copies_RB_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
# - blood 8 weeks highly correlated with LN 8 weeks
# - everything else not really highly correlated

## 3.3 Point IC50 and POS_AUC
pairs(dat_log %>% select(log_point_ic50_0_weekspost_ATI:pos_auc_8_weeks_post_ATI), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
# - a LOT of them are highly correlated
# - eps. measures in consecutive weeks
# - and IC50 & pos_auc in the same week

## 3.4 Viral loads and RNA copies
pairs(dat_log %>% select(log_peak_vl,log_vl_treat,log_RNA_copies_blood_8:log_RNA_copies_RB_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
# - VL at treatment highly correlated with RNA_blood_8, RNA_LN_8
# - RNA_blood_8 and RNA_LN_8 also highly correlated
# - peak VL highly correlated with RNA_blood_8

## 3.5 Viral loads and IC50
pairs(dat_log %>% select(log_peak_vl,log_vl_treat,log_point_ic50_0_weekspost_ATI:log_point_ic50_8_weekspost_ATI), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Not really correlated...

## (also checked antibody concentration and IC50)
pairs(dat_log %>% select(log_peak_gp41:log_gp120_treat,log_point_ic50_0_weekspost_ATI:log_point_ic50_8_weekspost_ATI), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Not really correlated either...
## some correlation between peak&treatment GP41 and IC50_8_weeks_post (.58 & .57)

## 3.6 RNA copies and IC50
pairs(dat_log %>% select(log_RNA_copies_blood_8:log_RNA_copies_RB_56,log_point_ic50_0_weekspost_ATI:log_point_ic50_8_weekspost_ATI), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Not really correlated...

# 4. Look at single-predictor models
Response = "Surv(rebound_time_days_post_ati, observed)"
All_covars = names(dat_log)[6:24]

# 4.1 check concordance, the c-statistic

C_stats = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  C_v = concordance(f, data=dat_log, timewt = "n")$concordance
  C_stats = c(C_stats, C_v)
}

C_stats = data.frame(Predictor = All_covars, Concordance = C_stats)

## show it with descending rank
C_stats %>% arrange(desc(Concordance))

# Predictor Concordance
# 1        pos_auc_0_weeks_post_ATI   0.7954545
# 2  log_point_ic50_0_weekspost_ATI   0.7727273
# 3  log_point_ic50_8_weekspost_ATI   0.7500000
# 4                  log_gp41_treat   0.7045455
# 5  log_point_ic50_2_weekspost_ATI   0.7045455
# 6  log_point_ic50_4_weekspost_ATI   0.7045455
# 7                   log_peak_gp41   0.6818182
# 8        pos_auc_4_weeks_post_ATI   0.6818182
# 9        pos_auc_2_weeks_post_ATI   0.6477273
# 10       pos_auc_8_weeks_post_ATI   0.6363636
# 11        log_RNA_copies_blood_56   0.5909091
# 12            log_RNA_copies_LN_8   0.5681818
# 13           log_RNA_copies_LN_56   0.5454545
# 14                 log_peak_gp120   0.5227273
# 15         log_RNA_copies_blood_8   0.5227273
# 16                log_gp120_treat   0.5113636
# 17                    log_peak_vl   0.4772727
# 18                   log_vl_treat   0.4772727
# 19           log_RNA_copies_RB_56   0.4545455


# 4.2 Cox Proportional Hazards model
## one example 
phmod = coxph(Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI,
              data = dat_log)
phmod_summ = summary(phmod)

