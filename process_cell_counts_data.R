# 06/15/2020
# read in and process cell counts data for B2 group
# Also: read in and process the "A01" status (MHC class labels)

# 06/21/2020: modify data error! RQc19 censor time 1/17/2020
# 98 days --> 168 days!

# 06/23/2020: read in the "animal_information" data


library(tidyverse)
setwd("~/Documents/Research_and_References/HIV_rebound_summer2020/")

# read in the CD4 and CD8 cell counts data
CD4 = read_csv("B1B2_Abs_CD4_flowCBC.csv")[,1:4]

CD8 = read_csv("B1B2_Abs_CD8_flowCBC.csv")[,1:4]

CD4$animal_id == CD8$animal_id # those labels are all the same...

cell_counts = cbind(CD4, CD8[,2:4])

# the data include B1 and B2, but only extract B2
dat_log = readRDS("reboundB2_logTrans_withVLAUC.rds")

B2_ids = dat_log$animal_id # all ordered alphabetically...

cell_ids = CD4$animal_id[which(CD4$animal_id %in% B2_ids)]

## then only need to order the data columns by the correct id order
cell_counts = cell_counts %>% filter(animal_id %in% B2_ids)
cell_counts = cell_counts[order(cell_ids),]

## check the order again
cell_counts$animal_id == B2_ids

# create log-transformed counts
cell_counts = cell_counts %>% 
  mutate_at(2:7, list(log=log10)) %>% 
  rename_at(vars( contains( "_log") ), list( ~paste("log", gsub("_log", "", .), sep = "_") ))

# save the cell counts data (including log transformed columns)
saveRDS(cell_counts,"B2_cell_counts.rds")

# combine the log counts with the big data table
log_cell_counts = cell_counts %>%
  select(log_Abs_CD4_week0:log_Abs_CD8_week8)

dat_log = cbind(dat_log, log_cell_counts)


# read in the A01 info
A01 = read_csv("animal_A01.csv")

## check animal IDs are all correct
A01$animal_id == B2_ids

# then combine with the big data table
dat_log = cbind(dat_log, A01$A01)
names(dat_log)[43] = "A01"

# save the big data table
saveRDS(dat_log, "reboundB2_logTrans_CellCounts_A01.rds")

# 06/21/2020: modify data error! RQc19 censor time 1/17/2020
# 98 days --> 168 days!

dat_log$rebound_time_days_post_ati[5] = 168
saveRDS(dat_log, "reboundB2_logTrans_CellCounts_A01.rds")

# 06/23/2020
# add animal info columns (Sex, challenge_time, dosage)
animal_info = read_csv("animal_information.csv")
animal_info$animal_id == dat_log$animal_id

# combine and save
dat_log = cbind(dat_log, animal_info[,3:5])
saveRDS(dat_log, "reboundB2_logTrans_CellCounts_AnimalInfo.rds")
