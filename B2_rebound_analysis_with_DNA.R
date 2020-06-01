# 06/01/2020
# B2 rebound survival analysis
# with DNA copies and viral_load_2

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

## 3.1: DNA vs RNA in blood
pairs(dat_log %>% 
        select(log_RNA_copies_blood_8,log_RNA_copies_blood_56, 
               log_DNA_copies_Blood_8:log_DNA_copies_Blood_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - RNA and DNA highly correlated on week 8 (also between weeks 8 and 16)
### - Not much correlation on week 56????

## 3.2: DNA vs RNA in lymph node
pairs(dat_log %>% 
        select(log_RNA_copies_LN_8,log_RNA_copies_LN_56, 
               log_DNA_copies_LN_8:log_DNA_copies_LN_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - RNA and DNA kinda correlated (0.55) on week 8
### - DNA copies highly correlated between weeks 8 and 16
### - But still, not much correlation for week 56 copies

## 3.3: DNA vs RNA in rectal biopsy
pairs(dat_log %>% 
        select(log_RNA_copies_RB_56, 
               log_DNA_copies_RB_16:log_DNA_copies_RB_56), 
      lower.panel = panel.cor, diag.panel = panel.hist)
### Findings:
### - Only see high correlation between RNA week 56 and DNA week 36
### - all the others are very low correlation
