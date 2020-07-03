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

# 1.2 fit a Cox PH model with `A01`
phmod_A01 = coxph(Surv(rebound_time_days_post_ati, observed)~A01, 
                   data = dat_log)
summary(phmod_A01)

### Summary:
# (Doesn't seem to be an important correlate at all!)
# Concordance= 0.477  (se = 0.11 )
# Likelihood ratio test= 0.08  on 1 df,   p=0.8
# Wald test            = 0.08  on 1 df,   p=0.8
# Score (logrank) test = 0.08  on 1 df,   p=0.8


# 2. check out the cell counts too
# 2.1 calculate concordance
# (re-do all the predictors to obtain the entire ranking)
Response = "Surv(rebound_time_days_post_ati, observed)"
All_covars = names(dat_log)[6:43]

C_stats = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  C_v = concordance(f, data=dat_log, timewt = "n")$concordance
  C_stats = c(C_stats, C_v)
}

C_stats = data.frame(Predictor = All_covars, Concordance = C_stats)

## show it with descending rank
C_stats %>% arrange(desc(Concordance))

### Not very important...
# log_Abs_CD4_week8 and log_Abs_CD8_week8 ranked at 12 and 13
# they are the only ones with >0.5 C-stats


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

## (using concordance, 5-fold)
cv_phmod_lasso2 = cv.glmnet(X, 
                            Surv(dat_log_sel$rebound_time_days_post_ati, 
                                 dat_log_sel$observed),
                            family = "cox", nfolds = 5,
                            type.measure = "C")
plot(cv_phmod_lasso2)
### also, somewhere between between 1 and 4 predictors: 
### C = close to 1

## extract coefficients using the 1st CV results
cv_phmod_lasso1$lambda.min
# [1] 0.3589552
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

# 4.3 Use Lasso results to obtain a "predictor inclusion ranking"
pred_in = NULL
coef_sign = NULL

for(l in cv_phmod_lasso1$lambda){
  coefficients = coef(phmod_lasso, s = l)
  active_index = which(coefficients != 0)
  active_coefficients = coefficients[active_index]
  active_predictors = attr(coefficients,"Dimnames")[[1]][active_index]
  
  #cat(active_predictors, "\n")
  
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

pred_in2 = NULL
LOO_deviance = NULL
for(l in cv_phmod_lasso3$lambda){
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


## put together a summary table of this thing
predictor_inclusion = data.frame(Inclusion_Rank = c(1:9),
                                 Predictor = pred_in,
                                 Coefficient = coef_sign,
                                 Rebound_Effect = rebound_effect,
                                 #Concordance = C_stats)
                                 Cum_concordance = Cum_C_stats,
                                 LOOCV_deviance = LOO_deviance)
predictor_inclusion
