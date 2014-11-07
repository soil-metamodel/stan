library(rstan);

Ctotal = 7.7;      # total carbon [mg C/g total]
gamma <- 0.9944;   # proportion of carbon in pool 1 at init
y0 <- Ctotal * c(gamma, 1 - gamma);  # initial pools

t0 <- 0.0;          # initial time
T <- 42;            # last time (days)
ts <- 1:T;          # solution times (1 day intervals)

k1 <- 0.1401;       # pool 1 decomp rate
k2 <- 0.3457;       # pool 2 decomp rate
alpha21 <- 0.8529;  # pool 2 to 1 transfer coeff
alpha12 <- 0.3742;  # pool 1 to 2 transfer coeff

data_vars <-  c("k1", "k2", "alpha21", "alpha12", 
                "y0", "t0",                       
                "T", "ts");                        

# deterministic simulation config
fit <- stan("soil_incubation_sim.stan",
            algorithm = "Fixed_param",
            data = data_vars,
            warmup=0, iter=1, chains=1)

fit_ss <- extract(fit);           
y_hat <- fit_ss$y_hat[1,,];           # simulated values
CO2_hat <- y_hat[,1] + y_hat[,2];     # total CO2
gamma_hat <- y_hat[,1] / CO2_hat;         # proportion CO2 in pool 1
eCO2_hat <- Ctotal - CO2_hat;

library('ggplot2');
df <- data.frame(t = c(0,ts), eCO2 = c(0,eCO2_hat));
evolved_co2_plot <- 
  ggplot(df, aes(x=t, y=eCO2)) + 
  geom_line() +
  scale_x_continuous(name="Days") +
  scale_y_continuous(name="Evolved CO2 (mgc g-1 soil)")
plot(evolved_co2_plot);
