library(rstan);
library(reshape2);
library(ggplot2);

totalC_t0 <- 7.7;         # total initial carbon (mg C g-1 C)

# parameter posterior mean estimates
# following Mueller & Sierra SoilR vignette (2014), p. 5

k1 <- 0.1401;             # pool 1 decomp rate
k2 <- 0.3457;             # pool 2 decomp rate
alpha21 <- 0.8529;        # pool 2 to 1 transfer coeff
alpha12 <- 0.3742;        # pool 1 to 2 transfer coeff
gamma <- 0.9345;          # partitioning coefficient
var_eCO2 <- 0.1256;       # observation noise variance

t0 <- 0.0;          # initial time
T <- 42;            # last time (days)
ts <- 1:T;          # solution times (1 day intervals)

data_vars <-  c("k1", "k2", "alpha21", "alpha12", "var_eCO2", 
                "totalC_t0", "gamma",
                "t0", "T", "ts");

# deterministic simulation config
fit <- stan("soil_incubation_sim.stan",
            algorithm = "Fixed_param",
            data = data_vars,
            warmup=0, iter=1, chains=1)
fit_ss <- extract(fit);  

eCO2_sim <- fit_ss$eCO2[1,];     # take first (and only) draw
eCO2_hat <- fit_ss$eCO2_hat[1,]; 

df <- data.frame(t = c(0,ts),
                 eCO2_hat = c(0, eCO2_hat),
                 eCO2 = c(0,eCO2_sim));

df.long <- melt(df, id.vars="t");

evolved_CO2_plot <-
  ggplot(df.long, aes(t,value,linetype=variable)) + 
  geom_line() + 
  scale_linetype(name="value",
                 labels=c("estimated","measured")) + 
  labs(x="days", 
       y="evolved CO2 (mgc g-1 soil)") +
  ggtitle("Evolved CO2 Simulation\n(2 pools, feedback)")
plot(evolved_CO2_plot);

