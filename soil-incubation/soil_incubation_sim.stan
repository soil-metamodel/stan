functions {
  real[] two_pool_feedback(real t, real[] y, real[] theta,
                           real[] x_r, int[] x_i) {
    real k1;
    real k2;
    real alpha21;
    real alpha12;

    real dy_dt[2];

    k1 <- theta[1];
    k2 <- theta[2];
    alpha21 <- theta[3];
    alpha12 <- theta[4];
    
    dy_dt[1] <- -k1 * y[1] + alpha12 * k2 * y[2];
    dy_dt[2] <- alpha21 * k1 * y[1] - k2 * y[2];

    return dy_dt;
  }
}
data {
  // model parameters
  real<lower=0> k1;       // pool 1 decomposition rate
  real<lower=0> k2;       // pool 2 decomposition rate
  real<lower=0> alpha21;  // transfer coeff from pool 2 to 1
  real<lower=0> alpha12;  // transfer coeff from pool 1 to 2

  // initial state and time
  real<lower=0> y0[2];
  real t0;

  // solution times
  int<lower=0> T;          // number of solution times
  real<lower=t0> ts[T];    // solution times
}
transformed data {
  // no real or integer data constants for the ODE
  real x_r[0];
  int x_i[0];
}
model {
}
generated quantities {
  real y_hat[T,2];   // solutions to initial value problem
  { 
    real theta[4];   // local var for parameter vector
    theta[1] <- k1;
    theta[2] <- k2;
    theta[3] <- alpha21;
    theta[4] <- alpha12;
    y_hat <- integrate_ode(two_pool_feedback,
                           y0, t0, ts, theta, x_r, x_i);
  }
                         
}
