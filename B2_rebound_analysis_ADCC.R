# 11/05/2020
# B2 rebound analysis with ADCC and ADCP data added
# as supplement to previous analysis

# 11/09/2020
# more stuff done:
# 1. analysis without that not-rebounded animal
# 2. analysis only focused on antibody responses (look at GP41 and GP120?)

library(tidyverse)
library(survival)
library(survminer)
library(glmnet)

setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

## latest data version
dat_log = readRDS('reboundB2_logTrans_ADCC_ADCP.rds')

## 11/09/20 tried: only focus on those who did rebound (sensitivity analysis)
#dat_log = dat_log %>% filter(observed == 1)


# 0. Visualize rebound time v.s. autogolous neutralization

dat_log$auto_neut = factor(c(0,0,1,0,1,1,0,1,1,0), levels=c(1,0),
                           labels=c('Yes','No'))

autoNeut_KM = survfit(Surv(rebound_time_days_post_ati, observed)~auto_neut, 
                     data=dat_log)
ggsurvplot(autoNeut_KM,  
           data = dat_log, 
           censor.shape="|", censor.size = 4,
           size = 1.5, #palette = c("#E7B800", "#2E9FDF"),
           xlab='Days after ATI', ylab='No-rebound probability',
           conf.int = FALSE,
           ggtheme = theme_bw(base_size=14))

# 1. Martingale test to check (non)linearity between 
#    survival and ADCC/ADCP
# 1.6: Use Martingale residuals to check (non)linearity
#      and select proper functional forms for the predictors
#      (only do this for RNA and IC50 covariates in the dataset)
Response = "Surv(rebound_time_days_post_ati, observed)"
Vars = names(dat_log)[47:48]
for(v in Vars){
  log_v = paste0("log(",v,")")
  log_v_2 = paste0("I(log(",v,")^2)")
  sqrt_v = paste0("sqrt(",v,")")
  f = as.formula(paste(paste(Response,v,sep = " ~ "),
                       log_v, log_v_2, sqrt_v, sep=" + "))
  print(ggcoxfunctional(f, data = dat_log))
}
## Notes on the "ggcoxfunctional" function:
## Plots Martingale residuals against transformed variables
## Should expect a linear trend of the points if linearity is satisfied
## Well, transformations don't seem to make a difference
## So right now, will just use the original scale
## But maybe create a log-scale ADCC??


# 2. check out "ADCC" at week 56 (just before ATI)
#   and ADCP at week 62 (at ATI) as well

## add a log_ADCC_week56 column
dat_log = dat_log %>% 
  mutate(log_ADCC_week56 = log(ADCC_week56))


# Univariate Concordance

Response = "Surv(rebound_time_days_post_ati, observed)"
All_covars = names(dat_log)[6:49]
## seems that log_ADCC_week56 doesn't make a difference by much
# so will keep on the original scale
All_covars = names(dat_log)[6:48]

C_stats = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  C_v = concordance(f, data=dat_log, timewt = "n")$concordance
  C_stats = c(C_stats, C_v)
}

C_stats = data.frame(Predictor = All_covars, Concordance = C_stats)

## show it with descending rank
C_stats %>% arrange(desc(Concordance))

## get a table of concordance
library(xtable)
C_stats_top = C_stats %>% arrange(desc(Concordance)) %>%
  slice(1:8)
xtable(C_stats_top, digits = 4)

## only look at ADCC and ADCP
C_stats %>% filter(Predictor %in%  c('ADCC_week56', 'ADCP_week62'))
# ADCC: 0.523 (ranked #25 among 43)
# ADCP: 0.500 (ranked #28 among 43)


## survplot for AUC0
phmod_auc0 = coxph(Surv(rebound_time_days_post_ati, observed) ~ 
                        pos_auc_0_weeks_post_ATI, 
                        data = dat_log)
summary(phmod_auc0)

## survival curves at representative values of AUC0
## vary AUC_0_week
auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean

auc0_fit = survfit(phmod_auc0, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Pos AUC \n0 weeks\npost ATI",
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           ggtheme = theme_bw(base_size = 14))


# 3. Try out "Coxnet": using the `glmnet` package
# 3.1 get the design matrix 
# (ignore A01 for now - can only deal with numeric values)
X = dat_log[,6:42] %>% as.matrix()

# 3.2 fit Lasso (alpha=1)
phmod_lasso = glmnet(X, 
                     Surv(dat_log$rebound_time_days_post_ati, dat_log$observed),
                     family = "cox")

# (cross validation)
## (using partial likelihood, 5-fold)
cv_phmod_lasso1 = cv.glmnet(X, 
                            Surv(dat_log$rebound_time_days_post_ati, dat_log$observed),
                            family = "cox", nfolds = 5)
plot(cv_phmod_lasso1)
### somewhere between 1 and 4 predictors...
### BUT! CV with 10 observations isn't very reliable
### Probably can only get "promising predictors"

## (using concordance, 5-fold)
cv_phmod_lasso2 = cv.glmnet(X, 
                            Surv(dat_log$rebound_time_days_post_ati, dat_log$observed),
                            family = "cox", nfolds = 5,
                            type.measure = "C")
plot(cv_phmod_lasso2)
### doesn't seem very helpful: using more predictors always yield higher C-index!

## extract coefficients using the 1st CV results
cv_phmod_lasso1$lambda.min
# [1] 0.215188
coefficients <- coef(phmod_lasso, s = cv_phmod_lasso1$lambda.min)
active_index <- which(coefficients != 0)
active_coefficients <-coefficients[active_index]
active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]

### put together a table for this result
cv_res = data.frame(Predictor = active_coefficients, 
                    Coefficient = active_coefficients)


# 3.3 Use Lasso results to obtain a "predictor inclusion ranking"
pred_in = NULL
coef_sign = NULL

for(l in cv_phmod_lasso1$lambda){
  coefficients = coef(phmod_lasso, s = l)
  active_index = which(coefficients != 0)
  active_coefficients = coefficients[active_index]
  active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]
  
  if(any(!active_predictors %in% pred_in)){
    new_index = which(!active_predictors %in% pred_in)
    pred_in = c(pred_in, active_predictors[new_index])
    new_coef_sign = ifelse(active_coefficients[new_index] > 0, 
                           "+", "-")
    coef_sign = c(coef_sign, new_coef_sign)
  }
}

## get a vector of "effect on rebound"
## (+: accelerate; -: delay)
rebound_effect = sapply(coef_sign, 
                        function(x) ifelse(x=="+","accelerate","delay")) %>%
  as.vector()

## put together a summary table of this thing
predictor_inclusion = data.frame(Inclusion_Rank = c(1:9),
                                 Predictor = pred_in,
                                 Coefficient = coef_sign,
                                 Rebound_Effect = rebound_effect)


# 4. Try "Lasso" Cox PH model with pre-selected predictors
# (to decrease colinearity)

# 4.1 put together a dataset with "manually selected" predictors
dat_log_sel = dat_log %>% 
  select(rebound_time_days_post_ati, observed,
         log_peak_vl_2, log_gp41_treat,
         log_RNA_copies_blood_56, log_RNA_copies_LN_56,
         log_RNA_copies_RB_56, log_DNA_copies_Blood_36,
         log_DNA_copies_LN_8, log_DNA_copies_LN_56,
         log_DNA_copies_RB_16, log_DNA_copies_RB_56,
         pos_auc_0_weeks_post_ATI,
         log_Abs_CD4_week0, log_Abs_CD4_week8,
         Challenge_times, Dosage)
X = dat_log_sel[,3:17] %>% as.matrix()
## here: 15 predictors

## (06/29/2020)
## a version with more predictors allowed
## cross corr < 0.8
dat_log_sel = dat_log %>% 
  select(rebound_time_days_post_ati, observed,
         log_peak_vl_2, log_vl_treat,
         log_peak_gp41, log_peak_gp120,
         log_RNA_copies_blood_56, log_RNA_copies_LN_56,
         log_RNA_copies_RB_56, 
         log_DNA_copies_Blood_16:log_DNA_copies_Blood_56,
         log_DNA_copies_LN_36, log_DNA_copies_LN_56,
         log_DNA_copies_RB_16:log_DNA_copies_RB_56,
         pos_auc_0_weeks_post_ATI, pos_auc_4_weeks_post_ATI,
         pos_auc_8_weeks_post_ATI,
         log_Abs_CD4_week0:log_Abs_CD4_week8,
         Challenge_times, Dosage)

## right now: 23 predictors
X = dat_log_sel[,3:25] %>% as.matrix()

## (11/5/2020)
# add in ADCC and ADCP!
dat_log_sel = dat_log %>% 
  select(rebound_time_days_post_ati, observed,
         log_peak_vl_2, log_vl_treat,
         log_peak_gp41, log_peak_gp120,
         log_RNA_copies_blood_56, log_RNA_copies_LN_56,
         log_RNA_copies_RB_56, 
         log_DNA_copies_Blood_16:log_DNA_copies_Blood_56,
         log_DNA_copies_LN_36, log_DNA_copies_LN_56,
         log_DNA_copies_RB_16:log_DNA_copies_RB_56,
         pos_auc_0_weeks_post_ATI, pos_auc_4_weeks_post_ATI,
         pos_auc_8_weeks_post_ATI,
         log_Abs_CD4_week0:log_Abs_CD4_week8,
         Challenge_times, Dosage, 
         ADCC_week56, ADCP_week62)

## now: 25 predictors in the pool
X = dat_log_sel[,3:27] %>% as.matrix()

# (11/09/2020)
# only include antibody response variables
# (also convert Sex and A01 to numeric variables)
dat_log_sel = dat_log %>% 
  mutate(sex = ifelse(Sex == 'F', 0, 1),
         MHC_A01 = ifelse(A01 == 'Pos', 1, 0)) %>%
  select(rebound_time_days_post_ati, observed, 
         Challenge_times, Dosage, 
         sex, MHC_A01, 
         log_peak_gp41, log_peak_gp120,
         log_gp41_treat, log_gp120_treat,
         log_point_ic50_0_weekspost_ATI:log_point_ic50_8_weekspost_ATI,
         pos_auc_0_weeks_post_ATI:pos_auc_8_weeks_post_ATI,
         ADCC_week56, ADCP_week62)

## here: 18 predictors in the pool
X = dat_log_sel[,3:20] %>% as.matrix()


# 4.2 fit Lasso (alpha=1)
phmod_lasso = glmnet(X, 
                     Surv(dat_log_sel$rebound_time_days_post_ati, 
                          dat_log_sel$observed),
                     family = "cox")

# (cross validation)
## (using partial likelihood, 5-fold)
cv_phmod_lasso1 = cv.glmnet(X, 
                            Surv(dat_log_sel$rebound_time_days_post_ati, 
                                 dat_log_sel$observed),
                            family = "cox", nfolds = 5)
plot(cv_phmod_lasso1)
### somewhere between 1 and 4 predictors...
### BUT! CV with 10 observations isn't very reliable
### Probably can only get "promising predictors"

### New on 11/05/2020: pretty much the same, 1~4 predictors

## (using concordance, 5-fold)
cv_phmod_lasso2 = cv.glmnet(X, 
                            Surv(dat_log_sel$rebound_time_days_post_ati, 
                                 dat_log_sel$observed),
                            family = "cox", nfolds = 5,
                            type.measure = "C")
plot(cv_phmod_lasso2)
### also, somewhere between between 1 and 4 predictors: 
### C = close to 1

### New on 11/05/2020: 2 predictors seem best

## extract coefficients using the 1st CV results
cv_phmod_lasso1$lambda.min
# [1] 0.3589552
# new 11/05/2020: 
# [1] 0.2254343
coefficients <- coef(phmod_lasso, s = cv_phmod_lasso1$lambda.min)
active_index <- which(coefficients != 0)
active_coefficients <-coefficients[active_index]
active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]

### put together a table for this result
cv_res = data.frame(Predictor = active_predictors, 
                    Coefficient = active_coefficients)
cv_res
# Predictor  Coefficient
# 1           log_gp41_treat -0.003705671
# 2 pos_auc_0_weeks_post_ATI -6.054194923

## New on 11/05/2020
#                 Predictor    Coefficient
# 1  log_RNA_copies_blood_56  -0.04274380
# 2 pos_auc_0_weeks_post_ATI -11.48880921
# 3        log_Abs_CD4_week8  -0.08234868
# 4          Challenge_times  -0.02683582
# 5              ADCC_week56  -0.06270668


# 4.3 Use Lasso results to obtain a "predictor inclusion ranking"
pred_in = NULL
coef_sign = NULL

for(l in phmod_lasso$lambda){
  coefficients = coef(phmod_lasso, s = l)
  active_index = which(coefficients != 0)
  active_coefficients = coefficients[active_index]
  active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]
  
  #cat(active_predictors, "\n")
  
  if(any(!active_predictors %in% pred_in)){
    cat('Added new predictors! Now active:', active_predictors, '\n')
    
    new_index = which(!active_predictors %in% pred_in)
    pred_in = c(pred_in, active_predictors[new_index])
    new_coef_sign = ifelse(active_coefficients[new_index] > 0, 
                           "+", "-")
    coef_sign = c(coef_sign, new_coef_sign)
  }
}

## 11/05/2020
## when lambda changes, we may end up selecting slightly different combos
## of variables than before (so it's not exactly cumulative)
## so only focus on all the 9 variables that eventually got selected with 
## min lambda

coef_sign = coef_sign[pred_in %in% active_predictors]
pred_in = pred_in[pred_in %in% active_predictors]

## get a vector of "effect on rebound"
## (+: accelerate; -: delay)
rebound_effect = sapply(coef_sign, 
                        function(x) ifelse(x=="+","accelerate","delay")) %>%
  as.vector()

## also get a vector of (univariate) concordance
Response = "Surv(rebound_time_days_post_ati, observed)"
All_covars = pred_in

C_stats = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  C_v = concordance(f, data=dat_log_sel, timewt = "n")$concordance
  C_stats = c(C_stats, C_v)
}

## 06/29/2020
## get a vector of multivariate (cumulative concordance)
Response = "Surv(rebound_time_days_post_ati, observed)"

Cum_C_stats = NULL
for(i in 1:length(pred_in)){
  covars = pred_in[1:i]
  f = as.formula(paste(Response,paste(covars,collapse = "+"),sep = " ~ "))
  mod = coxph(f, data=dat_log_sel)
  C_covars = mod$concordance['concordance'] %>% as.numeric()
  Cum_C_stats = c(Cum_C_stats, C_covars)
}

## 07/02/2020:
## get leave-one-out CV results
## (in terms of deviance, i.e., partial likelihood)
cv_phmod_lasso3 = cv.glmnet(X, 
                            Surv(dat_log_sel$rebound_time_days_post_ati, 
                                 dat_log_sel$observed),
                            family = "cox", nfolds = 10, keep = T)

## 11/09/2020:
## do stuff with only 9 data points (excluding RQc19)
cv_phmod_lasso3 = cv.glmnet(X, 
                            Surv(dat_log_sel$rebound_time_days_post_ati, 
                                 dat_log_sel$observed),
                            family = "cox", nfolds = 9, keep = T)

pred_in2 = NULL
LOO_deviance = NULL
for(l in phmod_lasso$lambda){
  coefficients = coef(phmod_lasso, s = l)
  active_index = which(coefficients != 0)
  #active_coefficients = coefficients[active_index]
  active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]
  
  #cat(active_predictors, "\n")
  
  if(any(!active_predictors %in% pred_in2)){
    new_index = which(!active_predictors %in% pred_in2)
    pred_in2 = c(pred_in2, active_predictors[new_index])
    
    l_index = which(cv_phmod_lasso3$lambda == l)
    LOO_deviance = c(LOO_deviance, cv_phmod_lasso3$cvm[l_index])
  }
}

## 11/05/2020
## same as before...
## only focus on all the 9 variables that eventually got selected with 
## min lambda
LOO_deviance = LOO_deviance[pred_in2 %in% active_predictors]
pred_in2 = pred_in2[pred_in2 %in% active_predictors]


## put together a summary table of this thing
## (use 1:8 for the 9 data point version)
predictor_inclusion = data.frame(#Inclusion_Rank = c(1:9),
                                 Inclusion_Rank = c(1:7),
                                 Predictor = pred_in,
                                 Coefficient = coef_sign,
                                 Rebound_Effect = rebound_effect,
                                 #Concordance = C_stats)
                                 Cum_concordance = Cum_C_stats,
                                 LOOCV_deviance = LOO_deviance)
predictor_inclusion

## get a latex version of this table
library(xtable)
## only include the columns that my slide can hold
predictor_inclusion  = predictor_inclusion %>%
  select(Predictor, 
         Effect = Rebound_Effect, 
         Cum.Concord = Cum_concordance,
         LOOCV_deviance)


xtable(predictor_inclusion, digits = 3)

## only get top 7 rows for 9 data point 
xtable(predictor_inclusion[1:7,], digits = 3)


## 07/13/2020
## add a bit more visualization on the "chosen" model out of Lasso + Cox PH
phmod_auc0_gp41 = coxph(Surv(rebound_time_days_post_ati, observed) ~ 
                          pos_auc_0_weeks_post_ATI + log_peak_gp41, 
                        data = dat_log)
summary(phmod_auc0_gp41)
### Summary:
# Concordance= 0.864  (se = 0.061 )
# Likelihood ratio test= 10.21  on 2 df,   p=0.006
# Wald test            = 6.62  on 2 df,   p=0.04
# Score (logrank) test = 8.74  on 2 df,   p=0.01

## 11/05/2020
## the "chosen" model now includes ADCC_week56
phmod_auc0_gp41_adcc = coxph(Surv(rebound_time_days_post_ati, observed) ~ 
                               pos_auc_0_weeks_post_ATI + log_peak_gp41 +
                               ADCC_week56, 
                             data = dat_log)
summary(phmod_auc0_gp41_adcc)

## (11/09/2020)
## plot coefficients and 95% CIs
coefs = coef(phmod_auc0_gp41_adcc)
names(coefs)[1] = 'pos_auc_0'
coef_dat = data.frame(variable=names(coefs), coefficient = coefs)
coef_dat = cbind(coef_dat, confint(phmod_auc0_gp41_adcc))
names(coef_dat)[3:4] = c('LB',"UB")

ggplot(data=coef_dat, aes(x=variable, y=coefficient)) +
  geom_hline(yintercept = 0, size=1.5, color='darkgray') +
  geom_errorbar(aes(ymin = LB, ymax = UB), size=1.2,
                color='darkblue', width=0.2) +
  geom_point(size=2, color='red') +
  labs(x='predictor') +
  coord_flip() +
  theme_bw(base_size = 14)

## try out only AUC_0 and ADCC_week56??
phmod_auc0_adcc = coxph(Surv(rebound_time_days_post_ati, observed) ~ 
                          pos_auc_0_weeks_post_ATI + ADCC_week56, 
                        data = dat_log)
summary(phmod_auc0_adcc)

# a bit more "significant" in terms of Wald test? 
## but ONLY AUC0 is significant at 0.05 level in either model...


## b) survival curves at representative values of each variable
## i) fix peak GP41 and ADCC_week56 at mean, vary AUC_0_week
mean(dat_log$log_peak_gp41)
# [1] 4.512378
mean(dat_log$ADCC_week56)
# [1] 23.73512
auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean

auc0_fit = survfit(phmod_auc0_gp41_adcc, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values,
                                        log_peak_gp41 = mean(dat_log$log_peak_gp41),
                                        ADCC_week56 = mean(dat_log$ADCC_week56)))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Pos AUC \n0 weeks\npost ATI",
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           caption = "Fix log_peak_gp41 at 4.5, ADCC_week56 at 23.7 (mean)",
           ggtheme = theme_bw(base_size = 14))

## ii) fix AUC_0_week and ADCC_week56 at mean, vary log peak gp41
sort(dat_log$log_peak_gp41)
# [1] 3.929253 4.147115 4.153393 4.159597 4.495433 4.612036 4.826567 4.836934
# [9] 4.846614 5.116840
mean(dat_log$pos_auc_0_weeks_post_ATI)
# [1] 0.24949
peakgp41_values = c(4, 4.5, 5) # 4.5 is approximately the mean value

peakgp41_fit = survfit(phmod_auc0_gp41_adcc, 
                       newdata = data.frame(pos_auc_0_weeks_post_ATI = mean(dat_log$pos_auc_0_weeks_post_ATI),
                                            ADCC_week56 = mean(dat_log$ADCC_week56),
                                            log_peak_gp41 = peakgp41_values))
ggsurvplot(peakgp41_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Peak GP41 \nconcentration\n(log-scale)",
           legend.labs=c("4.0", 
                         "4.5(mean)",
                         "5.0"),
           caption = "Fix pos_auc_0 at 0.25, ADCC_week56 at 23.7 (mean)",
           ggtheme = theme_bw(base_size = 14))

## iii) Fix AUC_0 and peak_gp41 at mean values, and vary ADCC_week56
sort(dat_log$ADCC_week56)
# [1] 18.65487 18.65487 21.36407 21.36407 21.54801 23.74600 24.35697 25.83522
# [9] 30.22471 31.60245
mean(dat_log$ADCC_week56)
ADCC_values = c(18.5, 23.7, 26, 31.5)

ADCC_fit = survfit(phmod_auc0_gp41_adcc, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = mean(dat_log$pos_auc_0_weeks_post_ATI),
                                        ADCC_week56 = ADCC_values,
                                        log_peak_gp41 = mean(dat_log$log_peak_gp41)))
ggsurvplot(ADCC_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "ADCC at \nWeek56\n(before ATI)",
           legend.labs=c("18.5", 
                         "23.7(mean)",
                         "26",
                         "31.5"),
           caption = "Fix pos_auc_0 at 0.25, log_peak_gp41 at 4.5 (mean)",
           ggtheme = theme_bw(base_size = 14))


## 11/09/2020
# the model with only antibody response variables
phmod_auc0_gp41treat = coxph(Surv(rebound_time_days_post_ati, observed) ~ 
                          pos_auc_0_weeks_post_ATI + log_gp41_treat, 
                        data = dat_log)
summary(phmod_auc0_gp41treat)

## plot coefficients and 95% CIs
coefs = coef(phmod_auc0_gp41treat)
names(coefs)[1] = 'pos_auc_0'
coef_dat = data.frame(variable=names(coefs), coefficient = coefs)
coef_dat = cbind(coef_dat, confint(phmod_auc0_gp41treat))
names(coef_dat)[3:4] = c('LB',"UB")

ggplot(data=coef_dat, aes(x=variable, y=coefficient)) +
  geom_hline(yintercept = 0, size=1.5, color='darkgray') +
  geom_errorbar(aes(ymin = LB, ymax = UB), size=1.2,
                color='darkblue', width=0.2) +
  geom_point(size=3, color='red') +
  labs(x='predictor') +
  coord_flip() +
  theme_bw(base_size = 14)

## visualization of survival curves
### i) fix log_gp41_treat
mean(dat_log$log_gp41_treat)
# [1] 4.491903

auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean

auc0_fit = survfit(phmod_auc0_gp41treat, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values,
                                        log_gp41_treat = mean(dat_log$log_gp41_treat)))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Pos AUC \n0 weeks\npost ATI",
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           caption = "Fix log_gp41_treat at 4.5(mean)",
           ggtheme = theme_bw(base_size = 14))

## ii) fix pos_auc_0_weeks_post_ATI
sort(dat_log$log_gp41_treat)

gp41_treat_values = c(3.9, 4.5, 5.1)

gp41_treat_fit = survfit(phmod_auc0_gp41treat, 
                         newdata = data.frame(pos_auc_0_weeks_post_ATI = mean(dat_log$pos_auc_0_weeks_post_ATI),
                                              log_gp41_treat = gp41_treat_values))

ggsurvplot(gp41_treat_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "GP41 (log-scale)\nconcentration\nat treatment",
           legend.labs=c("3.9", 
                         "4.5(mean)",
                         "5.1"),
           caption = "Fix pos_auc_0 at 0.25(mean)",
           ggtheme = theme_bw(base_size = 14))
