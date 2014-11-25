# script to fit 2-pool with feedback model to
# AK_T25 soil incubation data (included)

library(rstan);
library(ggplot2);


source("AK_T25.data.R");  

fit <- stan("soil_incubation.stan",
            data=c("totalC_t0", "t0", "N_t", "ts", "eCO2mean"),
            control=list(adapt_delta=0.95,
                         stepsize=0.01),
            chains=2, iter=100, seed=1234);

# to generate version of figure 2 in Mueller & Sierra (2014)
pairs_plot <- function(fit) {
  pairs(fit, pars=c("k1","k2","alpha21","alpha12","gamma"));
}

# to generate something like figure 4 in Mueller & Sierra (2014)
# technically:
#   * points are measurements
#   * line is posterior median for predicted evolved CO2
#   * interval is line +/- 1.96 posterior mean sd_eCO2 in the noise model
# you'd get slightly wider with proper posterior predictive
# using generated quantities to generate new data sets and taking 95%
# ints
pp_plot <- function(fit) {
  fit_ss <- extract(fit);
  sigma_hat <- median(fit_ss$sigma);
  eCO2_hat <- rep(NA,N_t);
  for (t in 1:N_t) {
    eCO2_hat[t] <- quantile(fit_ss$eCO2_hat[,t], c(0.5))
  }
  df_post <- data.frame(list(ts = ts, 
                             eCO2meas = eCO2mean,
                             eCO2_hat = eCO2_hat));
  ggplot(df_post, aes(x = ts)) +
    geom_ribbon(aes(ymin = eCO2_hat - 1.96 * sigma_hat,
                    ymax = eCO2_hat + 1.96 * sigma_hat),
                    fill="lightyellow") +
    geom_line(aes(y=eCO2_hat),colour="darkred") +
    geom_point(aes(y=eCO2meas),colour="darkblue") +
    labs(x="days", 
         y="evolved CO2 (mgC g-1 soil)") +
    ggtitle("Soil Incubation: 2 pools, feedback\n(posterior median, predictive 95%)");

   #  95% central interval for prediction
   #  geom_ribbon(aes(ymin=V3,ymax=V5),fill="lightyellow")
}

# this just shows what the sampler's up to over time
traces_plot <- function(fit) {
 traceplot(fit,c("k1","k2","alpha21","alpha12","gamma"))
}

