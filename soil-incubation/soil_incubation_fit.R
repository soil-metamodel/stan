# script to fit 2-pool with feedback model to
# AK_T25 soil incubation data (included)

library(rstan);

source("AK_T25.data.R");  

fit <- stan("soil_incubation.stan",
            data=c("totalC_t0", "t0", "T", "ts", "eCO2"),
            control=list(max_treedepth=7), 
            chains=1, iter=400, refresh=5);
