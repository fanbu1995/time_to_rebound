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

