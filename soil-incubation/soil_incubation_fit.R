# script to fit 2-pool with feedback model to
# AK_T25 soil incubation data (included)

library(rstan);
library(ggplot2);


source("AK_T25.data.R");  

fit <- stan("soil_incubation.stan",
            data=c("totalC_t0", "t0", "T", "ts", "eCO2"),
            control=list(adapt_delta=0.95,
                         stepsize=0.01),
            chains=1, iter=2000, refresh=1, seed=1234);

# to generate version of figure 2 in Mueller & Sierra (2014)
pairs_plot <- function(fit) {
  pairs(fit, pars=c("k1","k2","alpha21","alpha12","gamma"), condition="acc");
}

# to generate something like figure 4 in Mueller & Sierra (2014)
# technically:
#   * points are measurements
#   * line is posterior median for predicted evolved CO2
#   * interval is line +/- 1.96 posterior mean sd_eCO2 in the noise model
# you'd get slightly wider with proper posterior predictive
# using generated quantities to generate new data sets and taking 95% ints
pp_plot <- function(fit) {
  fit_ss <- extract(fit);
  sd_eCO2_hat <- median(fit_ss$sd_eCO2);
  eCO2_hat <- matrix(NA,T,3);
  for (t in 1:T) {
    eCO2_hat[t,] <- quantile(fit_ss$eCO2_hat[,t], c(0.025,0.5,0.975));
  }
  df <- data.frame(cbind(ts,eCO2,eCO2_hat));

  ggplot(df, aes(x=ts)) +
    geom_ribbon(aes(ymin = V4 - 1.96 * sd_eCO2_hat, ymax = V4 + 1.96 * sd_eCO2_hat), fill="lightyellow") +
    geom_line(aes(y=V4),colour="darkred") +
    geom_point(aes(y=eCO2),colour="darkblue") +
    labs(x="days", 
         y="evolved CO2 (mgc g-1 soil)") +
    ggtitle("Soil Incubation: 2 pools, feedback\n(posterior median, predictive 95%)");

   #  95% central interval for prediction
   #  geom_ribbon(aes(ymin=V3,ymax=V5),fill="lightyellow")
}

# this just shows what the sampler's up to over time
traces_plot <- function(fit) {
 traceplot(fit,c("k1","k2","alpha21","alpha12","gamma"))
}

