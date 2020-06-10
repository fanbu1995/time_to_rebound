# 06/05/2020
# look a bit more at post ATI viral levels
# and see what factors influence them

# 06/09/2020
# extract post-rebound:
# - peak viral load
# - ared-under-VL-curve


# load required packages and set working directory
library(tidyverse)
library(ggplot2)
library(zoo)
library(ggpubr)


setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")


VL = readRDS("B2_viral_load.rds")

# visualize by animal IDs
ggplot(data=VL, aes(x=week_infection, y=log_viral_load)) + 
  geom_vline(xintercept = c(8,62), size=1) +
  geom_line(aes(color=animal_id)) +
  labs(x="Weeks post infection", 
       y="Log viral load",
       color="Animal ID")+
  scale_x_continuous(breaks = c(0,8,25,50,62,75)) +
  theme_bw(base_size = 14)

# visualize only the part after ATI
VL_post = VL %>% filter(days_post_ati >= 0)

# load rebound time data as well (mark up rebound time in plots)
dat_log = readRDS("reboundB2_logTrans_withVLAUC.rds")
rebound_times = dat_log %>% 
  select(animal_id, rebound_time_days_post_ati)
rebound_times$rebound_time_days_post_ati[rebound_times$animal_id=="RQc19"] = NA

ggplot(data=VL_post, aes(x=days_post_ati, y=log_viral_load)) +
  geom_hline(yintercept = log(60, base=10), 
             size = 0.5, linetype = "dashed") +
  geom_vline(aes(xintercept = rebound_time_days_post_ati),
             data = rebound_times, 
             size = 0.5, linetype = "dashed") +
  geom_line(aes(color=animal_id)) +
  geom_point(aes(color=animal_id), size=1) +
  labs(x="Days post ATI", 
       y="Viral load (log10)",
       color="Animal ID")+
  #scale_x_continuous(breaks = c(0,8,25,50,62,75)) +
  scale_y_continuous(limits = c(0,8)) +
  theme_bw(base_size = 14)+
  facet_wrap(~animal_id, ncol = 5)

## (function for AUC)
# a function to calculate AUC first
get_AUC <- function(vl, times){
  times = times[-which(is.na(vl))]
  vl = vl[-which(is.na(vl))]
  id = order(times)
  sum(diff(times[id])*zoo::rollmean(vl[id],2))
}

# extract post rebound viral info
VL_post_rebound = VL_post %>% 
  filter(animal_id != "RQc19")

## impute values below cutoff as **half** of the viral load
VL_post_rebound$viral_load[which(VL_post_rebound$viral_load <= 60)] = 30
VL_post_rebound$log_viral_load[which(VL_post_rebound$viral_load <= 60)] = log(30, base=10)

VL_post_rebound = VL_post_rebound %>%
  group_by(animal_id) %>%
  summarise(log_peakVL_postATI = max(log_viral_load, na.rm = T),
            log_VL_AUC_postATI = get_AUC(log_viral_load, week_infection))

## combine it with `dat_log` and save it
dat_log_rebound = dat_log %>% filter(animal_id != "RQc19")
dat_log_rebound = cbind(dat_log_rebound, VL_post_rebound[,2:3])

saveRDS(dat_log_rebound, "rebound_B2_with_postATI_VL.rds")

# 1) look at peak vs AUC
ggplot(data = dat_log_rebound,aes(y=log_VL_AUC_postATI, x=log_peakVL_postATI)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_bw(base_size = 14)

## plot with correlation and p-values
ggscatter(data = dat_log_rebound, 
          x = "log_peakVL_postATI", y = "log_VL_AUC_postATI",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = T,
          cor.coef = T)

# 2) look at peak/AUC vs all predictors
## compile a table of correlation coefficients
corrs_peak = NULL
corrs_AUC = NULL
for(n in All_covars){
  corrs_peak = c(corrs_peak, 
                 cor(dat_log_rebound$log_peakVL_postATI,
                     dat_log_rebound[,n]))
  corrs_AUC = c(corrs_AUC, 
                 cor(dat_log_rebound$log_VL_AUC_postATI,
                     dat_log_rebound[,n]))
}

sort(corrs_peak)
sort(corrs_AUC)

### A little bit of negative correlation, but mostly positive correlation

All_covars[order(corrs_peak, decreasing = T)]
# top 6
# [1] "log_vl_treat"                   "log_DNA_copies_Blood_8"        
# [3] "log_RNA_copies_LN_8"            "log_RNA_copies_blood_8"        
# [5] "log_vl_auc"                     "log_DNA_copies_Blood_36" 

All_covars[order(corrs_AUC, decreasing = T)]
# top 6
# [1] "log_vl_treat"                   "log_RNA_copies_LN_8"           
# [3] "log_DNA_copies_Blood_36"        "log_DNA_copies_Blood_8"        
# [5] "log_DNA_copies_Blood_16"        "log_RNA_copies_blood_8"        


### Important factors: 1) pre-ART viral load, 2) CA-DNA&RNA

# 3) visualize a bit
## post peak vs VL at treatment
ggscatter(data = dat_log_rebound, 
          y = "log_peakVL_postATI", x = "log_vl_treat",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = T,
          cor.coef = T)

### (the outlier: RJm19 -- low postATI VL, but relatively high pre-ART VL)

## post peak vs log_DNA_copies_Blood_8 (DNA in blood, at treat)
ggscatter(data = dat_log_rebound, 
          y = "log_peakVL_postATI", x = "log_DNA_copies_Blood_8",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = T,
          cor.coef = T)

### (but here RJm10 not really an outlier)

## post VL AUC vs log_vl_treat
ggscatter(data = dat_log_rebound, 
          y = "log_VL_AUC_postATI", x = "log_vl_treat",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = T,
          cor.coef = T)

## post VL AUC vs log_RNA_copies_LN_8
ggscatter(data = dat_log_rebound, 
          y = "log_VL_AUC_postATI", x = "log_RNA_copies_LN_8",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = T,
          cor.coef = T)

### (outlier: RLg19 overly high post VL AUC)


