functions {

  /**
   * ODE system for two pool model with feedback and no inputs.
   *
   * System State C is two dimensional with C[1] and C[2] 
   * being carbon in pools 1 and 2.
   *
   * The system has parameters
   *
   *   theta = (k1, k2, alpha21, alpha12)
   *
   * where
   *
   *   k1:       pool 1 decomposition rate
   *   k2:       pool 2 decomposition rate
   *   alpha21:  transfer coefficient from pool 2 to pool 1
   *   alpha12:  transfer coefficient from pool 1 to pool 2
   *
   * The system time derivatives are
   *
   *   d.C[1] / d.t  =  -k1 * C[1]  +  alpha12 * k2 * C[2]
   *
   *   d.C[2] / d.t  =  alpha21 * k1 * C[1]  -  k2 * C[2]
   *
   * @param t time at which derivatives are evaluated.
   * @param C system state at which derivatives are evaluated.
   * @param theta parameters for system.
   * @param x_r real constants for system (empty).
   * @param x_i integer constants for system (empty).
   */
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
    
    dC_dt[1] <- -k1 * C[1] + alpha12 * k2 * C[2];
    dC_dt[2] <- alpha21 * k1 * C[1] - k2 * C[2];

    return dC_dt;
  }
}
data {
  real<lower=0> totalC_t0;   // initial total carbon

  real t0;                   // initial time
  int<lower=0> T;            // number of measurement times
  real<lower=t0> ts[T];      // measurement times

  real<lower=0> eCO2[T];     // measured cumulative evolved carbon
}
transformed data {
  real x_r[0];               // no real data for ODE system
  int x_i[0];                // no integer data for ODE system
}
parameters {
  real<lower=0> k1;          // pool 1 decomposition rate
  real<lower=0> k2;          // pool 2 decomposition rate
  real<lower=0> alpha21;     // transfer coeff from pool 2 to 1
  real<lower=0> alpha12;     // transfer coeff from pool 1 to 2
  real<lower=0,upper=1> gamma;  // partitioning coefficient
  real<lower=0> sd_eCO2;     // observation std dev
}
model {
  real C0[2];                 // initial state
  real theta[4];              // ODE parameters
  real C_hat[T,2];            // predicted pool content

  C0[1] <- gamma * totalC_t0;
  C0[2] <- (1 - gamma) * totalC_t0;

  theta[1] <- k1;
  theta[2] <- k2;
  theta[3] <- alpha21;
  theta[4] <- alpha12;

  C_hat <- integrate_ode(two_pool_feedback, 
                         C0, t0, ts, theta, x_r, x_i);
  for (t in 1:T)
    eCO2[t] ~ normal(totalC_t0 - sum(C_hat[t]), sd_eCO2);

  k1 ~ cauchy(0,1);
  k2 ~ cauchy(0,1);
  alpha21 ~ cauchy(0,1);
  alpha12 ~ cauchy(0,1);
  sd_eCO2 ~ cauchy(0,1);
  gamma ~ beta(10,1);
}
