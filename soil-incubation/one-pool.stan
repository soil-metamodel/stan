functions {

  real K(real k0, real Q10, real T) {
    return k0 * pow(Q10, (T - 15)/10);
  }

  real[] model_one(real t, real[] Cs, real[] theta,
                   real[] x_r, int[] x_i) {
    real k0;
    real Q10;

    real T;

    k0 <- theta[1];
    Q10 <- theta[2];

    T <- x_r[1];

    real[1] dCs_dt;
    dCs_dt[1] <- K(k0,Q10,T) * Cs[1];
    return dCs_dt;
  }
}
data {
  real<lower=0> Cs_t0;
  int<lower=0> N;         // num observation times
  real<lower=0> t[N];     // observation times
  real<lower=0> temp[N];  // temperatutes
  real<lower=0> Cs[T];    // soil C observations
}
transformed data {
  real[1] x_r;
  int[0] x_i;
  
  x_r[1] <- 
}

parameters {
  real<lower=0> Q10;
  real<lower=0> k0;
}
model {
  k0 ~ normal(0,10);
  Q10 ~ normal(0,10);

  {
    real[T,1] Cs_hat;
    Cs_hat <- integrate_ode(model_one, 
                            Cs_t0, t0, ts, theta, x_r, x_i);
  }  
}
