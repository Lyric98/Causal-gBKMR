library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

## File name: generate_popn_3t_highconf_linearLY_binary_probit.R

## Load necessary packages
## (If you are missing them, install them: install.packages("MASS"))
library(MASS)

## Set seed for reproducibility
set.seed(123)

## 1) Simulate a large "population" size from which you can later sample
N <- 30000

## 2) Baseline confounders: e.g., sex (binary) and baseline waist0
sex <- rbinom(N, size = 1, prob = 0.5)      # 0/1
waist0 <- 90 + 15 * rnorm(N)                # baseline waist circumference

## 3) Baseline exposures: logM1_0, logM2_0, logM3_0 depend on baseline confounders
##    (sex, waist0) plus random error. You can tune the coefficients below.
##    For correlation among the baseline exposures, you can draw the errors
##    from a multivariate normal. For simplicity, we show independent normals.
err0 <- mvrnorm(n = N, mu = c(0,0,0),
                Sigma = diag(3))  # 3x3 identity => no correlation
e1_0 <- err0[,1]
e2_0 <- err0[,2]
e3_0 <- err0[,3]

logM1_0 <-  0.50 * sex + 0.01 * waist0 + e1_0
logM2_0 <-  0.30 * sex + 0.02 * waist0 + e2_0
logM3_0 <- -0.10 * sex + 0.04 * waist0 + e3_0

## 4) Time-1 waist: waist1
##    Suppose waist1 depends on baseline waist0, baseline exposures, plus noise:
waist1 <- 20 +
  0.3 * waist0 +
  0.5 * sex +
  1.0 * (logM1_0 + logM2_0 + logM3_0) +    # depends on baseline exposures
  rnorm(N, sd = 1)

## 5) Time-1 exposures: logM1_1, logM2_1, logM3_1
##    Depend on (waist1 + baseline exposures). E.g.:
err1 <- mvrnorm(n = N, mu = c(0,0,0), Sigma = diag(3))
e1_1 <- err1[,1]
e2_1 <- err1[,2]
e3_1 <- err1[,3]

logM1_1 <-  0.2 * sex + 0.2 * waist1 +
  0.4 * (logM1_0 + logM2_0 + logM3_0) + e1_1
logM2_1 <-  0.1 * sex + 0.2 * waist1 +
  0.3 * (logM1_0 + logM2_0 + logM3_0) + e2_1
logM3_1 <- -0.1 * sex + 0.2 * waist1 +
  0.2 * (logM1_0 + logM2_0 + logM3_0) + e3_1

## 6) Time-2 waist: waist2
##    Suppose waist2 depends on waist1, baseline waist0, time-1 exposures, etc.
waist2 <- 30 +
  0.2 * waist0 +
  0.3 * waist1 +
  0.5 * sex +
  0.5 * (logM1_1 + logM2_1 + logM3_1) +
  rnorm(N, sd = 1)

## 7) Time-2 exposures: logM1_2, logM2_2, logM3_2
##    Depend on waist2 + waist1 + the time-1 exposures
err2 <- mvrnorm(n = N, mu = c(0,0,0), Sigma = diag(3))
e1_2 <- err2[,1]
e2_2 <- err2[,2]
e3_2 <- err2[,3]

logM1_2 <-  0.2 * sex + 0.3 * waist2 + 0.1 * waist1 +
  0.3 * (logM1_1 + logM2_1 + logM3_1) + e1_2
logM2_2 <-  0.1 * sex + 0.3 * waist2 + 0.1 * waist1 +
  0.2 * (logM1_1 + logM2_1 + logM3_1) + e2_2
logM3_2 <- -0.1 * sex + 0.3 * waist2 + 0.1 * waist1 +
  0.2 * (logM1_1 + logM2_1 + logM3_1) + e3_2

## 8) Generate a binary outcome Y under a probit model:
##    Y* = linear predictor + error, error ~ Normal(0,1); Y=1{Y*>0}.
##    Adjust coefficients below as you like (stronger/weaker effects).
lp <- -1 +
  0.3*sex +
  0.05*waist0 + 0.10*waist1 + 0.20*waist2 +
  0.05*(logM1_0 + logM2_0 + logM3_0) +
  0.05*(logM1_1 + logM2_1 + logM3_1) +
  0.05*(logM1_2 + logM2_2 + logM3_2)

ystar <- lp + rnorm(N, 0, 1)
Y     <- ifelse(ystar > 0, 1, 0)

## 9) Combine all data into one data frame; create ID
id <- seq_len(N)
sim_popn <- data.frame(
  id,
  sex,
  waist0, waist1, waist2,
  logM1_0, logM2_0, logM3_0,
  logM1_1, logM2_1, logM3_1,
  logM1_2, logM2_2, logM3_2,
  Y
)

## 10) Save to disk as RDS
saveRDS(sim_popn, file = "popn_3t_highconf_linearLY_binary_probit.rds")

message("Finished generating popn_3t_highconf_linearLY_binary_probit.rds")
