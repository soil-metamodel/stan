# script to extract AK_T25 soil incubation data from SoilR

library(SoilR);
library(rstan);

totalC_t0 <- 7.7;       # not included in data, so hard code here
t0 <- 0;
N_t <- 25;                # calculated by inspection of eCO2
ts <- eCO2[1:25,2];
eCO2mean <- eCO2[1:25,3];   
eCO2sd <- eCO2[1:25,4];     

stan_rdump(c("t0","N_t","ts","eCO2mean","eCO2sd","totalC_t0"),
           "AK_T25.data.R");
