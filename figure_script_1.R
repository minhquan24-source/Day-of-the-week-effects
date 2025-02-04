


# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)
library(RColorBrewer)

origin_date <- as.Date("2024-09-12")
date_column <- "notification_date"

# Loading required functions
source('R/ps_single_analysis_scripts.R')


######################################################################################
## Loading in the data 

df_cov <- read.csv(paste("data/SARSCOV2-case-count-", origin_date, ".csv", sep=""))


######################################################################################
## Formatting the data a little

df_cov <- df_cov[df_cov$test_type=="PCR",]
df_cov[,date_column] <- as.Date(df_cov[,date_column])
######################################################################################
##  Set limits on dates to consider

# Covid date limits
max_date <- origin_date
min_date <- as.Date("2023-01-01")
df_cov <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>=min_date,]

######################################################################################
## Formatting the dates a little


df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1

df_cov <- df_cov[order(df_cov$time_index),]

################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models

ps_single_mod <- stan_model('stan/ps_single_final.stan')



############################################################################################
## Reading in the stan models

cov_fit  <- readRDS(paste('fitted_stan_models/', 'cov_fit-overall.rds', sep=""))


###############################################################################################
## Extract smoothed estimates for COVID
cov_mod_inc_dow <- ps_single_incidence_dow(cov_fit,
                                           week_effect = 7,
                                           DOW = (df_cov$time_index %% 7)+1,
                                           X=df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])

cov_mod_inc <- ps_single_incidence(cov_fit, df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])
cov_mod_gr <- ps_single_growth_rate(cov_fit, df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])



#########################################################################################################
# Plotting figures
#########################################################################################################
first_date <- as.Date("2024-09-10") - 27*7
###########################################################
## Figure 1
###########################################################


cov1 <- ggplot(cov_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=df_cov, aes(x=notification_date, y=cases) , linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  #coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")


cov2 <- ggplot(cov_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  #coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(-0.052,0.052))+
  scale_y_continuous(breaks=c(-0.04,0.0,0.04))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")




cov1 + cov2 

ggsave('figure/Figure1.png', width=8, height=8)


