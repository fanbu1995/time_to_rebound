# 03/18/2021

# add GP41 and GP120 concentration at ATI timepoint

library(tidyverse)
library(readxl)

setwd('~/Documents/Research_and_References/HIV_rebound_summer2020/')

# load dataset (after adding week 8 measurements of ADCC and ADCP)
dat = readRDS('reboundB2_logTrans_032021.rds')

# read GP41 and GP120 data
elisa = read_excel('B2_ELISA_CH505gp120_MNgp41_remaining ATI.xlsx')
names(elisa)[4] = 'antigen'

elisa = elisa %>% 
  select(animal_id, antigen, week_infection, days_post_ATI,
         conc = conc_binding_abs)

# select ATI timepoints for GP41 and GP120
GP41 = elisa %>% filter(antigen == "HIV-1_MN gp41_Immuno_Dx",
                        days_post_ATI == 0) %>%
  select(animal_id, conc) %>%
  arrange(animal_id)

GP41$animal_id

## manually add in RVh19 and RWc19 data
RVh19_data = elisa %>% 
  filter(antigen == "HIV-1_MN gp41_Immuno_Dx", animal_id=='RVh19',
         days_post_ATI == 4) %>%
  select(conc) %>% pull()
## (this should be close enough to 0 days past ATI; 
##  don't have ATI measurement!!)

GP41_ext = data.frame(animal_id=c('RWc19', 'RVh19'), 
                      conc=c(1141.535, RVh19_data))
## NOTE: 1141.535 is from the GP41 excel file in Box

GP41 = rbind(GP41, GP41_ext) %>% arrange(animal_id)

## check animal_id
GP41$animal_id == dat$animal_id



GP120 = elisa %>% filter(antigen == "HIV CH505TF_08gp120/293F",
                         days_post_ATI == 0) %>%
  select(animal_id, conc)

# then manually add in data for RVh19 
#  (only have 4 days past ATI measurement)
GP120_ext = elisa %>% 
  filter(antigen == "HIV CH505TF_08gp120/293F", 
         animal_id=='RVh19', days_post_ATI==4) %>% 
  select(animal_id, conc)

GP120 = rbind(GP120, GP120_ext) %>% arrange(animal_id)

## check animal_id
GP120$animal_id == dat$animal_id


# combine with dataset
dat$log_gp41_ATI = log(GP41$conc, base = 10)
dat$log_gp120_ATI = log(GP120$conc, base = 10)

# save data
saveRDS(dat, 'reboundB2_logTrans_032021.rds')
