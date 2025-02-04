setwd("C:/Users/EALESO/R Projects/Winter 2024 paper")
# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)

origin_date <- as.Date("2024-09-12")
date_column <- "notification_date"

# Loading required functions
source('ps_single_analysis_scripts.R')

# Load the data 
df_cov <- read.csv(paste("data/SARSCOV2-case-count-", origin_date, ".csv", sep=""))
df_cov <- df_cov[df_cov$test_type=="PCR",]

# Set limits on dates to consider
max_date <- origin_date
min_date <- as.Date("2023-01-01")
df_cov <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>=min_date,]

df_cov[,date_column] <- as.Date(df_cov[,date_column])
df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1

df_cov <- df_cov[order(df_cov$time_index),]

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models

ps_single_mod <- stan_model('stan/ps_single_final.stan')

#############################################################################################################################################
## Fitting to COVID data (overall)

# Calculate the locations of equally spaced knots
knots <- get_knots(df_cov$time_index, days_per_knot = 5, spline_degree = 3)

cov_data <- list(num_data = nrow(df_cov),
                 num_knots = length(knots),
                 knots = knots,
                 spline_degree=3,
                 Y = df_cov$cases,
                 X = df_cov$time_index,
                 week_effect = 7,
                 DOW = (df_cov$time_index %% 7)+1 ) 

cov_fit <- sampling(ps_single_mod,
                    iter= 5000,
                    warmup = 1000,
                    chains=4,
                    data = cov_data)

saveRDS(cov_fit, paste('fitted_stan_models/', 'cov_fit-overall.rds', sep=""))


#############################################################################################################################################
## Fitting to COVID data (real-time analysis)

max_dates_considered <- max(df_cov$notification_date) - seq(0, 7*25, by=7)

for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1

  knots <- get_knots(df_tmp$time_index, days_per_knot = 5, spline_degree = 3)
  
  tmp_data <- list(num_data = nrow(df_tmp),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df_tmp$cases,
                   X = df_tmp$time_index,
                   week_effect = 7,
                   DOW = (df_tmp$time_index %% 7)+1 ) 
  
  tmp_fit <- sampling(ps_single_mod,
                      iter= 5000,
                      warmup = 1000,
                      chains=4,
                      data = tmp_data)
  
  saveRDS(tmp_fit, paste('fitted_stan_models/',max_date, '-cov_fit.rds', sep=""))
  
}