# script to extract AK_T25 soil incubation data from SoilR

library(SoilR);
library(rstan);

totalC_t0 <- 7.7;       # not included in data, so hard code here
t0 <- 0;
T <- 25;                # calculated by inspection of eCO2
ts <- eCO2[1:25,2];
eCO2 <- eCO2[1:25,3];   # nasty masking of SoilR's eCO2 variable

stan_rdump(c("t0","T","ts","eCO2","totalC_t0"),
           "AK_T25.data.R");
