# 06/01/2020
# B2 rebound survival analysis
# with DNA copies and viral_load_2

# 06/03/2020
# re-run some of the code with VL AUC included

# 06/18/2020
# re-do the univariate Cox PH model part
# and adjust for FDR (the B-H method)

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

## also read in DNA copies data
dat_DNA = read_csv("B2_CA_DNA.csv")
## log transform the DNA copies
dat_DNA = dat_DNA %>% 
  mutate_at(vars(contains("DNA")), list(log=log)) %>% 
  rename_at(vars( contains( "_log") ), list( ~paste("log", gsub("_log", "", .), sep = "_") ))
dat_DNA_log = dat_DNA %>% select(log_DNA_copies_Blood_8:log_DNA_copies_RB_56)
## check that the animal IDs match up
dat_DNA$animal_id == dat_log$animal_id
## then combine the log DNA copies with all other columns
dat_log = cbind(dat_log, dat_DNA_log)

## (save it again)
saveRDS(dat_log, "reboundB2_logTrans_withDNA.rds")

# 2. check correlation of DNA counts against RNA counts
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

## 2.1: DNA vs RNA in blood
pairs(dat_log %>% 
        select(log_RNA_copies_blood_8,log_RNA_copies_blood_56, 
               log_DNA_copies_Blood_8:log_DNA_copies_Blood_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - RNA and DNA highly correlated on week 8 (also between weeks 8 and 16)
### - Not much correlation on week 56????

## 2.2: DNA vs RNA in lymph node
pairs(dat_log %>% 
        select(log_RNA_copies_LN_8,log_RNA_copies_LN_56, 
               log_DNA_copies_LN_8:log_DNA_copies_LN_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - RNA and DNA kinda correlated (0.55) on week 8
### - DNA copies highly correlated between weeks 8 and 16
### - But still, not much correlation for week 56 copies

## 2.3: DNA vs RNA in rectal biopsy
pairs(dat_log %>% 
        select(log_RNA_copies_RB_56, 
               log_DNA_copies_RB_16:log_DNA_copies_RB_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - Only see high correlation between RNA week 56 and DNA week 36
### - all the others are very low correlation

## 2.4: also check viral loads vs. DNA
pairs(dat_log %>% select(log_peak_vl, log_vl_treat, 
                         log_DNA_copies_Blood_8:log_DNA_copies_RB_56),
      lower.panel = panel.cor, diag.panel = panel.hist)


# 06/03/2020
# (load the latest data version)
dat_log = readRDS("reboundB2_logTrans_withVLAUC.rds")


# 3. Univariate survival models
Response = "Surv(rebound_time_days_post_ati, observed)"
All_covars = names(dat_log)[6:35]

# (06/03/2020)
All_covars = names(dat_log)[6:36]

# (06/18/2020)
All_covars = names(dat_log)[6:42]

# 3.1 check concordance, the c-statistic

C_stats = NULL

for(v in All_covars){
  f = as.formula(paste(Response,v,sep = " ~ "))
  C_v = concordance(f, data=dat_log, timewt = "n")$concordance
  C_stats = c(C_stats, C_v)
}

C_stats = data.frame(Predictor = All_covars, Concordance = C_stats)

## show it with descending rank
C_stats %>% arrange(desc(Concordance))

### Findings:
### - Pretty much the same as before
### - the DNA counts don't have high concordance
### - the new "peak viral load" measure has even lower concordance

### (06/03/2020)
### - Still, pretty much same as before
### - log_vl_auc does better than other VL measures, ranked at #14

# 3.2 Univariate Cox Proportional hazards model

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

all_res_dat$predictor = All_covars
all_res_dat

## (06/18/2020)
## adjust by FDR (q-values proposed by B and H)
LR_pvals = all_res_dat$lik_ratio_test
LR_qvals = p.adjust(LR_pvals, method = "fdr")

## pick out the ones with FDR <= 0.2
below_thres = which(LR_qvals <= 0.2)
All_covars[below_thres]
# [1] "pos_auc_0_weeks_post_ATI"
LR_qvals[below_thres]
# [1] 0.08808077



## save this result too
write.csv(all_res_dat,"time_to_rebound/univariate_CoxPH_results_all.csv",
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

### EXACTLY the same as before! (the same on 06/03/2020)
# predictor lik_ratio_test
# 1 log_point_ic50_0_weekspost_ATI    0.015338725
# 2 log_point_ic50_8_weekspost_ATI    0.017391100
# 3       pos_auc_0_weeks_post_ATI    0.002380561


## 06/08/2020
## Adjust for multiple testing (using Sime's inequality)
## Test if any p_(i) < i*0.05/n

nvar = length(All_covars)
Is = seq_along(All_covars)

lrt_pvals = all_res_dat$lik_ratio_test
Thres = 0.05/nvar * Is

sort(lrt_pvals) < Thres

### All of them are FALSE! So at least using this method, can't reject the null...


# 4. Bi-variate (two-predictor) Cox PH model

## Since pos_auc_0_weeks_post_ATI is highly correlated with 
## - all the other AUCs, and
## - all the pointi_ic50s
## we don't consider those predictors as the second one
## (i.e., only consider VLs, Antibodies, RNA copies, DNA copies)

## (update on 06/03/2020: 
## correlation with log_vl_auc = .457, not very high)

f_auc0 = "Surv(rebound_time_days_post_ati, observed) ~ pos_auc_0_weeks_post_ATI"

Vars = names(dat_log)[c(6:16,25:35)]

# 06/03/2020
Vars = names(dat_log)[c(6:16,25:36)]
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

add_pred_res %>% arrange(AIC_linear) %>% head()

### There is a change: incluing log_peak_vl_2 has an edge!
# second_predictor AIC_linear AIC_interation
# 1           log_peak_vl_2   18.01492       19.91185
# 2             log_peak_vl   19.07615       20.25586
# 3    log_RNA_copies_RB_56   21.85564       21.55774

### 06/03/2020: log_peak_vl_2 still the best, 
### but log_vl_auc not bad either
# second_predictor AIC_linear AIC_interation
# 1           log_peak_vl_2   18.01492       19.91185
# 2             log_peak_vl   19.07615       20.25586
# 3              log_vl_auc   21.34036       23.25653

## look at this new bivariate model
phmod_auc0_peakVL2 = coxph(update(as.formula(f_auc0), ~ . + log_peak_vl_2), 
                          data = dat_log)
summary(phmod_auc0_peakVL2)
### Summary:
### Concordance: 0.932 (same as before with log_peak_vl)
### pos_auc_0_weeks_post_ATI: negative effect on hazard (delays rebound)
### log_peak_vl: positive effect on hazard (accelerates rebound)
### ALTHOUGH none of the effects is significantly non-zero (same as before)

confint(phmod_auc0_peakVL2)
#                               2.5 %   97.5 %
# pos_auc_0_weeks_post_ATI -166.849287 12.46821
# log_peak_vl_2              -1.521619 13.41001


## a) baseline survival curve (at the mean values)
ggsurvplot(survfit(phmod_auc0_peakVL2, data=dat_log), 
           palette = c("#2E9FDF"),
           ggtheme = theme_bw())
# HUGE confidence intervals

## b) survival curves at representative values of each variable
## i) fix peak VL at mean, vary AUC_0_week
mean(dat_log$log_peak_vl_2)
# [1] 5.8631
auc0_values = c(0.15, 0.25, 0.3, 0.4) # 0.25 is approximately the mean

auc0_fit = survfit(phmod_auc0_peakVL2, 
                   newdata = data.frame(pos_auc_0_weeks_post_ATI = auc0_values,
                                        log_peak_vl_2 = mean(dat_log$log_peak_vl_2)))
ggsurvplot(auc0_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "POS AUC \n0 weeks\npost ATI",
           legend.labs=c("0.15", 
                         "0.25(mean)",
                         "0.30",
                         "0.40"),
           caption = "Fix log_peak_VL_2 at 5.86 (mean)",
           ggtheme = theme_bw(base_size = 14))

## ii) fix AUC_0_week at mean, vary log peak VL 2
sort(dat_log$log_peak_vl_2)
# [1] 4.366423 5.196136 5.597164 5.651816 5.812189 6.008954 6.103290
# [8] 6.469425 6.483587 6.942020
mean(dat_log$pos_auc_0_weeks_post_ATI)
# [1] 0.24949
peakVL_values = c(4.5,5.2,5.8,6.5) # 5.8 is approximately the mean value

peakVL_fit = survfit(phmod_auc0_peakVL2, 
                     newdata = data.frame(pos_auc_0_weeks_post_ATI = mean(dat_log$pos_auc_0_weeks_post_ATI),
                                          log_peak_vl_2 = peakVL_values))
ggsurvplot(peakVL_fit, conf.int = FALSE, 
           data = dat_log,
           legend = "right",
           legend.title = "Peak viral load 2\n(log-scale)",
           legend.labs=c("4.5", 
                         "5.2",
                         "5.8(mean)",
                         "6.5"),
           caption = "Fix pos_auc_0_weeks_post_ATI at 0.25 (mean)",
           ggtheme = theme_bw(base_size = 14))

## c) check proportional hazards assumption
ggcoxzph(cox.zph(phmod_auc0_peakVL2))
# global p-value = 0.0576
# AUC0 p-value = 0.363
# log peak VL p-value = 0.9432
# Not very strong evidence to refute proportional hazards assumption

## d) look at influential observations/outliers
##    (using deviance residuals)
ggcoxdiagnostics(phmod_auc0_peakVL2, type = "deviance",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
## deviance quite small --> no obvious outliers
