# 05/28/2020
# first look at the B2 rebound data (cleaned by Yishu Gong)
# try to reproduce previous analysis (and do a bit more)

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


# 1.6: Use Martingale residuals to check (non)linearity
#      and select proper functional forms for the predictors
#      (only do this for RNA and IC50 covariates in the dataset)
Response = "Surv(rebound_time_days_post_ati, observed)"
Vars = names(dat)[18:26]
for(v in Vars){
  log_v = paste0("log(",v,")")
  log_v_2 = paste0("I(log(",v,")^2)")
  sqrt_v = paste0("sqrt(",v,")")
  f = as.formula(paste(paste(Response,v,sep = " ~ "),
                       log_v, log_v_2, sqrt_v, sep=" + "))
  print(ggcoxfunctional(f, data = dat))
}
## Notes on the "ggcoxfunctional" function:
## Plots Martingale residuals against transformed variables
## Should expect a linear trend of the points if linearity is satisfied

## Findings:
## Checked log(), ^2, sqrt()
## There exists some non-linearity, and transformation doesn't 
## magically solve it, but log() does shrink the large range
## and generally makes the trend seem more linear


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

# (06/23/2020)
# add "log_peak_vl_2"
pairs(dat_log %>% select(log_peak_vl:log_gp120_treat, log_peak_vl_2), 
      lower.panel = panel.cor, diag.panel = panel.hist)

### log_peak_vl, log_peak_vl_2: cor = 0.99

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


# 4.2 Cox Proportional Hazards model (univariate)
## one example 
phmod = coxph(Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI,
              data = dat_log)
phmod_summ = summary(phmod)

## Do:
# For each variable, extract:
# coef, exp(coef)
# concordance, LRT p-val, Wald test p-val, Score (log-rank) test p-val
# AIC, BIC (smaller values --> better)

## Also (TO DO): some kind of small-sample adjustment on those hypothesis tests??

all_res = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  phmod = coxph(f, data=dat_log)
  phmod_summ = summary(phmod)
  this_v = c(phmod_summ$coefficients[,c("coef","exp(coef)")],
             phmod_summ$concordance[1],
             phmod_summ$logtest["pvalue"] %>% as.numeric(),
             phmod_summ$waldtest["pvalue"] %>% as.numeric(),
             phmod_summ$sctest["pvalue"] %>% as.numeric(),
             AIC(phmod), BIC(phmod))
  all_res = rbind(all_res, this_v)
}

all_res_dat = as.data.frame(all_res, row.names = All_covars)
names(all_res_dat) = c("coef","exp(coef)","concordance",
                       "lik_ratio_test", "wald_test",
                       "log_rank_test", "AIC", "BIC")
all_res_dat

all_res_dat$predictor = All_covars

## save the results in the local directory for later usage
write.csv(all_res_dat,"time_to_rebound/univariate_CoxPH_results.csv",
          row.names = FALSE)

## look at pvalue < .05 ones
## log_rank_test
all_res_dat %>% filter(log_rank_test < 0.05)
## Wald test
all_res_dat %>% filter(wald_test < 0.05)
## likelihood ratio test
## (said to have best behavior with small sample size)
all_res_dat %>% filter(lik_ratio_test < 0.05) %>% 
  select(predictor, lik_ratio_test)
# 1 log_point_ic50_0_weekspost_ATI
# 2 log_point_ic50_8_weekspost_ATI
# 3       pos_auc_0_weeks_post_ATI (this has the strongest evidence of non-zero effects)

## rank the single predictors by AIC
all_res_dat %>% arrange(AIC) %>% 
  select(predictor, AIC)

## Examine two univariate models
## 1) predictor: pos_auc_0_weeks_post_ATI
phmod_auc0 = coxph(Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI,
                   data = dat_log)

### 1.a: baseline survival curve (at the MEAN value of pos_auc_0_weeks_post_ATI)
### (should be the same as the plain KM curve...?)
ggsurvplot(survfit(phmod_auc0, data=dat_log), palette = c("#2E9FDF"),
           ggtheme = theme_bw())

### 1.b: survival curves at certain representative values of pos_auc_0_weeks_post_ATI
sort(dat_log$pos_auc_0_weeks_post_ATI)
# [1] 0.1097 0.1474 0.1511 0.2184 0.2512 0.2660 0.2887 0.2959 0.3293 0.4372
mean(dat_log$pos_auc_0_weeks_post_ATI)
# [1] 0.24949
auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean
#auc0_values = quantile(dat_log$pos_auc_0_weeks_post_ATI, c(0.25,0.5,0.75))

auc0_fit = survfit(phmod_auc0, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "POS AUC \n0 weeks\npost ATI",
           # legend.labs=c("AUC_0_weeks=0.15", 
           #               "AUC_0_weeks=0.25(mean)",
           #               "AUC_0_weeks=0.30",
           #               "AUC_0_weeks=0.40"),
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           ggtheme = theme_bw(base_size = 14))
# the confidence bands are HUGE though...

### 1.c: test on Schoenfeld residuals 
###      (to validate proportional hazard assumption)
ggcoxzph(cox.zph(phmod_auc0))

### There isn't a clear pattern of residuals accross time
### p-value = 0.6772 -> no strong evidence against the PH assumption

### Notes on this:
# The function cox.zph() provides a convenient solution to test the proportional hazards assumption for each covariate included in a Cox regression model fit.
# 
# For each covariate, the function cox.zph() correlates the corresponding set of scaled Schoenfeld residuals with time, to test for independence between residuals and time. Additionally, it performs a global test for the model as a whole.

### 1.d: look at influential observations/outliers
### (using deviance residuals)
### Note on this:
###  - Positive values: died too soon
###  - Negative values: lived too long
ggcoxdiagnostics(phmod_auc0, type = "deviance",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
### deviance not too large, though hard to tell (too few points)


## 2) predictor: log_point_ic50_8_weekspost_ATI
## (not using "log_point_ic50_0_weekspost_ATI" 
##  because it is highly correlated with "pos_auc_0_weeks_post_ATI")
phmod_ic50_8 = coxph(Surv(rebound_time_days_post_ati, observed) ~ log_point_ic50_8_weekspost_ATI,
                     data = dat_log)

### 2.a: baseline survival curve (at the MEAN value of log_point_ic50_8_weekspost_ATI)
### (still the same as the plain KM curve...?)
ggsurvplot(survfit(phmod_ic50_8, data=dat_log), 
           palette = c("#2E9FDF"),
           ggtheme = theme_bw())

### 2.b: survival curves at certain representative values of log_point_ic50_8_weekspost_ATI
sort(dat_log$log_point_ic50_8_weekspost_ATI)
# [1] 2.998244 3.064655 3.080653 3.087281 3.103271 3.292609 3.380890 3.402849
# [9] 3.941457 4.055358
mean(dat_log$log_point_ic50_8_weekspost_ATI)
# [1] 3.340727
point_ic50_8_values = c(3.00, 3.34, 3.7, 4.0) # 3.34 is approximately the mean

ic50_8_fit = survfit(phmod_ic50_8, 
                     newdata = data.frame(log_point_ic50_8_weekspost_ATI =
                                            point_ic50_8_values))
ggsurvplot(ic50_8_fit, conf.int = TRUE, 
           data = dat_log,
           legend = "right",
           legend.title = "Point IC50 \n8 weeks\npost ATI\n(log-scale)",
           legend.labs=c("3.00", 
                         "3.34(mean)",
                         "3.70",
                         "4.00"),
           ggtheme = theme_bw(base_size = 14))
# the confidence bands are HUGE though...

### 2.c: test on Schoenfeld residuals 
###      (to validate proportional hazard assumption)
ggcoxzph(cox.zph(phmod_ic50_8))

### There isn't a clear pattern of residuals accross time
### p-value = 0.8128 -> no strong evidence against the PH assumption

### 2.d: look at influential observations/outliers
### (using deviance residuals)
ggcoxdiagnostics(phmod_ic50_8, type = "deviance",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
### deviance not too large, though hard to tell (too few points)
### Observation 7 has deviance close to 2... (died a bit too soon)
### Animal RTp19, rebound time = 7 days


# 4.3 Cox Proportional Hazards model (bivariate)
#     (extend 4.2 with ONE additional predictor)

# 4.3.a: Add another predictor to "pos_auc_0_weeks_post_ATI"
#        (+ a linear term AND an interaction term)
#        (select models by AIC)

## Since pos_auc_0_weeks_post_ATI is highly correlated with 
## - all the other AUCs, and
## - all the pointi_ic50s
## we don't consider those predictors as the second one
## (i.e., only consider VLs, Antibodies, RNA copies)

f_auc0 = "Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI"

Vars = names(dat_log)[6:16]
AIC_linear = NULL
AIC_inter = NULL
for(v in Vars){
  v_inter = paste0("pos_auc_0_weeks_post_ATI * ",v)
  f_linear = as.formula(paste(f_auc0,v,sep = " + "))
  f_inter = as.formula(paste(f_auc0,v,v_inter,sep = " + "))
  phmod_v_linear = coxph(f_linear, data = dat_log)
  phmod_v_inter = coxph(f_inter, data = dat_log)
  AIC_linear = c(AIC_linear, AIC(phmod_v_linear))
  AIC_inter = c(AIC_inter, AIC(phmod_v_inter))
}

add_pred_res = data.frame(second_predictor = Vars, 
                          AIC_linear = AIC_linear,
                          AIC_interation = AIC_inter)

## using "log_peak_vl" as second predictor is the best choice by far
## linear term AIC = 19.076
## (with AUC_0_week only: AIC = 22.97863)

## look a bit more closely at the bivariate model
phmod_auc0_peakVL = coxph(update(as.formula(f_auc0), ~ . + log_peak_vl), 
                          data = dat_log)
summary(phmod_auc0_peakVL)
### Summary:
### Concordance: 0.932
### pos_auc_0_weeks_post_ATI: negative effect on hazard (delays rebound)
### log_peak_vl: positive effect on hazard (accelerates rebound)
### ALTHOUGH none of the effects is significantly non-zero

## a) baseline survival curve (at the mean values)
ggsurvplot(survfit(phmod_auc0_peakVL, data=dat_log), 
           palette = c("#2E9FDF"),
           ggtheme = theme_bw())
# HUGE confidence intervals

## b) survival curves at representative values of each variable
## i) fix peak VL at mean, vary AUC_0_week
mean(dat_log$log_peak_vl)
# [1] 6.005144
auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean

auc0_fit = survfit(phmod_auc0_peakVL, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values,
                                        log_peak_vl = mean(dat_log$log_peak_vl)))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "POS AUC \n0 weeks\npost ATI",
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           caption = "Fix log_peak_VL at 6.0 (mean)",
           ggtheme = theme_bw(base_size = 14))

## ii) fix AUC_0_week at mean, vary log peak VL
sort(dat_log$log_peak_vl)
# [1] 4.409933 5.225751 5.788897 5.872156 6.085302 6.158818 6.260071 6.515874
# [9] 6.522444 7.212188
mean(dat_log$pos_auc_0_weeks_post_ATI)
# [1] 0.24949
peakVL_values = c(5,5.5,6.0,7) # 6 is approximately the mean value

peakVL_fit = survfit(phmod_auc0_peakVL, 
                     newdata = data.frame(pos_auc_0_weeks_post_ATI = mean(dat_log$pos_auc_0_weeks_post_ATI),
                                          log_peak_vl = peakVL_values))
ggsurvplot(peakVL_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Peak viral load\n(log-scale)",
           legend.labs=c("5.00", 
                         "5.50",
                         "6.00(mean)",
                         "7.00"),
           caption = "Fix pos_auc_0_weeks_post_ATI at 0.25 (mean)",
           ggtheme = theme_bw(base_size = 14))

## c) check proportional hazards assumption
ggcoxzph(cox.zph(phmod_auc0_peakVL))
# global p-value = 0.06567
# AUC0 p-value = 0.1651
# log peak VL p-value = 0.9997
# Not very strong evidence to refute proportional hazards assumption

## d) look at influential observations/outliers
##    (using deviance residuals)
ggcoxdiagnostics(phmod_auc0_peakVL, type = "deviance",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
## deviance quite small --> no obvious outliers

# 4.3.b: exhaustive search of all bivariate models
#  (still: linear term AND interaction term)
Response = "Surv(rebound_time_days_post_ati, observed)"
n_vars = length(All_covars)

Predictor1 = NULL
Predictor2 = NULL
AIC_linear = NULL
AIC_inter = NULL
for(i in 1:(n_vars-1)){
  for(j in 2:n_vars){
    v1 = All_covars[i]
    v2 = All_covars[j]
    v_inter = paste0(v1," * ",v2)
    f_linear = as.formula(paste(Response, 
                                paste(v1,v2,sep = " + "),
                                sep = " ~ "))
    f_inter = as.formula(paste(Response, 
                               paste(v1,v2,v_inter,sep = " + "),
                               sep = " ~ "))
    phmod_v_linear = coxph(f_linear, data = dat_log)
    phmod_v_inter = coxph(f_inter, data = dat_log)
    AIC_linear = c(AIC_linear, AIC(phmod_v_linear))
    AIC_inter = c(AIC_inter, AIC(phmod_v_inter))
    Predictor1 = c(Predictor1, v1)
    Predictor2 = c(Predictor2, v2)
  }
}

two_pred_res = data.frame(Predictor1 = Predictor1,
                          Predictor2 = Predictor2,
                          AIC_linear = AIC_linear,
                          AIC_interation = AIC_inter)
head(two_pred_res)
## sort by AIC_linear
two_pred_res %>% arrange(AIC_linear) %>% head()

## sort by AIC_interaction
two_pred_res %>% arrange(AIC_interation) %>% head()

## SAME AS BEFORE:
## Best model is
## Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI + log_peak_vl
