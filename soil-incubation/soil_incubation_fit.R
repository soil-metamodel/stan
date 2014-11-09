# script to fit 2-pool with feedback model to
# AK_T25 soil incubation data (included)

library(rstan);

source("AK_T25.data.R");  

fit <- stan("soil_incubation.stan",
            data=c("totalC_t0", "t0", "T", "ts", "eCO2"),
# control=list(adapt_delta=0.9,
#                         stepsize=0.005),
#                         # metric="dense_e",
#                         # max_treedepth=9),
            chains=1, iter=4000, refresh=1, seed=1234);

# to generate figure 2 in Mueller & Sierra (2014)
# pairs(fit,condition="acc")

# to generate equivalent of figure 3 in Mueller & Sierra (2014)
