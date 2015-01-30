functions {
  real[] two_pool_feedback(real t, real[] C, real[] theta,
                           real[] x_r, int[] x_i) {
    real k1;
    real k2;
    real alpha21;
    real alpha12;
    real dC_dt[2];
    k1 <- theta[1];
    k2 <- theta[2];
    alpha21 <- theta[3];
    alpha12 <- theta[4];
    dC_dt[1] <- -k1 * C[1] + alpha21 * C[2];
    dC_dt[2] <- - k2 * C[2] + alpha12 * C[1] ;
    return dC_dt;
  }
  real[] evolved_CO2(int N_t, real t0, real[] ts,
                     real gamma, real totalC_t0,
                     real k1, real k2, 
                     real alpha21, real alpha12,
                     real[] x_r, int[] x_i) {
    real C_t0[2];
    real theta[4];
    real C_hat[N_t,2];
    real eCO2_hat[N_t];
    C_t0[1] <- gamma * totalC_t0;
    C_t0[2] <- (1 - gamma) * totalC_t0;
    theta[1] <- k1;
    theta[2] <- k2;
    theta[3] <- alpha21;
    theta[4] <- alpha12;
    C_hat <- integrate_ode(two_pool_feedback, 
                           C_t0, t0, ts, theta, x_r, x_i);
    for (t in 1:N_t)
      eCO2_hat[t] <- totalC_t0 - sum(C_hat[t]);
    return eCO2_hat;
  }
}
data {
  real<lower=0> totalC_t0;
  real t0;
  int<lower=0> N_t;
  real<lower=t0> ts[N_t];
  real<lower=0> eCO2mean[N_t];
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  real<lower=0> k1;
  real<lower=0> k2;
  real<lower=0> alpha21;
  real<lower=0> alpha12;
  real<lower=0,upper=1> gamma;
  real<lower=0> sigma;
}
transformed parameters {
  real eCO2_hat[N_t];
  eCO2_hat <- evolved_CO2(N_t, t0, ts, gamma, totalC_t0,
                          k1, k2, alpha21, alpha12, x_r, x_i);
}
model {
  gamma ~ beta(10,1);
  k1 ~ normal(0,1);
  k2 ~ normal(0,1);
  alpha21 ~ normal(0,1);
  alpha12 ~ normal(0,1);
  sigma ~ cauchy(0,1);
  eCO2mean ~ normal(eCO2_hat, sigma);
}
