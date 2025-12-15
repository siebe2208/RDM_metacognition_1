library(tidyverse)
library(cmdstanr)
library(posterior)
library(here)

##################################################################
# Data Cleaning
##################################################################

# Extract data
data = read.csv(here("Analysis", "RDM_reportz_sub0_2.csv"))

# Filter for main and different scale
data_main = data[data$Trialtype == "Main",]
data_main = data_main[data_main$scale == "conf", ]

# Transform variables for model
data_main$SR_conf = as.numeric(data_main$SR_conf)
data_main$SR_conf = (data_main$SR_conf + 1)/2

data_main <- data_main %>%
  mutate(coherence = ifelse(cor_resp == "down", -(coherence), coherence))
data_main$up = ifelse(data_main$resp == "up", 1, 0)

#Delete trials where NA
data_main = data_main[!is.na(data_main$SR_conf), ]
sum(is.na(data_main$SR_conf))
data_main = data_main[!is.na(data_main$resp), ]
sum(is.na(data_main$resp))

rows <- which(data_main$SR_conf < 0.2 & data_main$cor == 1)
data_main = data_main[-rows,]
##################################################################
# Extract data for model
##################################################################
N = nrow(data_main)
binom_y = data_main$cor
RT = data_main$RTdec
conf = data_main$SR_conf
minRT = min(data_main$RTdec)
x = data_main$coherence
ACC = data_main$cor

##################################################################
# Model Fit
##################################################################

data_list = list(N = N, RT = RT, Conf = conf, X = x, minRT = minRT, ACC = ACC, binom_y = binom_y)

model = cmdstan_model(here("Analysis", "single_sub_nocopula.stan"))

fit = model$sample(data = data_list, parallel_chains = 4)

fit$save_object(here("Analysis", "fit_pilot_siebe.rds"))
fit <- readRDS(here("Analysis", "fit_pilot_siebe.rds"))
