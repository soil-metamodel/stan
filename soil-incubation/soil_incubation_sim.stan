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
  // model parameters
  real<lower=0> k1;       // pool 1 decomposition rate
  real<lower=0> k2;       // pool 2 decomposition rate
  real<lower=0> alpha21;  // transfer coeff from pool 2 to 1
  real<lower=0> alpha12;  // transfer coeff from pool 1 to 2
  real<lower=0> var_eCO2; // observation square error

  real totalC_t0;         // initial total carbon
  real<lower=0,upper=1> gamma;  // partitioning coefficient

  real t0;                // initial time
  int<lower=0> T;         // number of solution times
  real<lower=t0> ts[T];   // solution times
}
transformed data {
  real x_r[0];            // no real data for ODE system
  int x_i[0];             // no integer data for ODE system
}
model {
  /* intentionally empty for simulation */
}
generated quantities {
  real C_hat[T,2];   // solutions to initial value problem
  real eCO2_hat[T];  // cumulative total evolved carbon
  real eCO2[T];      // cumulative total evolved carbon w. noise
  { 
    real theta[4];   // local var for parameters
    real C0[2];
    C0[1] <- gamma * totalC_t0;
    C0[2] <- (1 - gamma) * totalC_t0;
    theta[1] <- k1;
    theta[2] <- k2;
    theta[3] <- alpha21;
    theta[4] <- alpha12;
    C_hat <- integrate_ode(two_pool_feedback,
                           C0, t0, ts, theta, x_r, x_i);
    for (t in 1:T) {
      eCO2_hat[t] <- totalC_t0 - sum(C_hat[t]);
      // fmax(0,...) because content can't be negative
      eCO2[t] <- fmax(0,normal_rng(eCO2_hat[t], sqrt(var_eCO2)));
    }
  }
}
