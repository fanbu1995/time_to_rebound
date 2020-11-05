# 11/04/2020
# process ADCP and ADCC data

library(tidyverse)
library(readxl)

setwd('~/Documents/Research_and_References/HIV_rebound_summer2020/')

# load previous data as reference
dat = readRDS('reboundB2_logTrans_CellCounts_AnimalInfo.rds')

# load data
adcc = read_excel('20201020_ADCC_B2_tms.sjb_amara.fouda.xls')
adcp = read_excel('20201020_ADCP_B2_sjb_fouda.xls')

# 1. ADCC
adcc = adcc %>% filter(week_challenge == 56) %>%
  select(animal_id, 11) %>%
  arrange(animal_id)

# check animal ids match up
adcc$animal_id == dat$animal_id

# change col name of ADCC
names(adcc)[2] = 'ADCC_week56'
glimpse(adcc)

# 2. ADCP
adcp = adcp %>% filter(week_ati == 0) %>%
  select(animal_id, 11) %>%
  arrange(animal_id)

# check animal ids match up
adcp$animal_id == dat$animal_id

# change col name of ADCP
names(adcp)[2] = 'ADCP_week62'
glimpse(adcp)


# 3. attach those new columns to the previous dataset
dat$ADCC_week56 = adcc$ADCC_week56
dat$ADCP_week62 = adcp$ADCP_week62

# 4. save it 
saveRDS(dat, 'reboundB2_logTrans_ADCC_ADCP.rds')
