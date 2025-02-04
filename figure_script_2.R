


# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)

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

######################################################################################
##  Set limits on dates to consider

# Covid date limits
max_date <- origin_date
min_date <- as.Date("2023-01-01")
df_cov <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>=min_date,]


######################################################################################
## Formatting the dates a little
df_cov[,date_column] <- as.Date(df_cov[,date_column])

df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1


df_cov <- df_cov[order(df_cov$time_index),]

################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Loading Stan models

ps_single_mod <- stan_model('stan/ps_single_final.stan')


############################################################################################
## Reading in the stan models and getting modelled estimates
############################################################################################
max_dates_considered <- max(df_rsv$notification_date) - seq(0, 7*25, by=7)

#############################################
## COVID Stan models
for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  cov_fit <- readRDS(paste('fitted_stan_models/',max_date, '-cov_fit.rds', sep=""))
  
  
  tmp_mod_inc <- ps_single_incidence(cov_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  tmp_mod_gr <- ps_single_growth_rate(cov_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  
  tmp_mod_inc$max_date <- max_date
  tmp_mod_gr$max_date <- max_date
  
  tmp_mod_inc$col <- i%%2
  tmp_mod_gr$col <- i%%2
  
  if(i==1){
    cov_mod_inc <- tmp_mod_inc
    cov_mod_gr <- tmp_mod_gr
    
  } else{
    cov_mod_inc <- rbind(cov_mod_inc, tmp_mod_inc)
    cov_mod_gr <- rbind(cov_mod_gr, tmp_mod_gr)
  }
  
}



##########################################################################################################

##################################################
# COVID
cov_mod_inc$min_date <- cov_mod_inc$max_date-7
cov_mod_inc$mask <- 0.0
cov_mod_inc[cov_mod_inc$time>=cov_mod_inc$min_date &cov_mod_inc$time<=cov_mod_inc$max_date,]$mask <- 1.0

cov_mod_gr$min_date <- cov_mod_gr$max_date-7
cov_mod_gr$mask <- 0.0
cov_mod_gr[cov_mod_gr$time>=cov_mod_gr$min_date &cov_mod_gr$time<=cov_mod_gr$max_date,]$mask <- 1.0

cov_mod_gr$col <- as.factor(cov_mod_gr$col)
cov_mod_inc$col <- as.factor(cov_mod_inc$col)


############################################################################################
## Making the figures
############################################################################################

############################################################################################
# Figure 4
############################################################################################
first_date <- as.Date("2024-09-10") - 25*7
cols <- c("blue4","red4")

############################################
## COVID panels
############################################

cov1 <- ggplot(cov_mod_inc[cov_mod_inc$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, color=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=cov_mod_inc[cov_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=cov_mod_inc[cov_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")


cov2<- ggplot(cov_mod_gr[cov_mod_gr$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, color=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  geom_line(data=cov_mod_gr[cov_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=cov_mod_gr[cov_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  theme_bw(base_size = 14)+
  ylab("Growth rate")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")


cov1 + cov2  + plot_layout(nrow=2, heights=c(2,1))

ggsave('figure/Figure4.png', width=8, height=10)


