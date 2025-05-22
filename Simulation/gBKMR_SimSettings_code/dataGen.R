library(lubridate)
library(dplyr)
library(data.table)
library(ggplot2) 
library(haven)
library(readstata13) 
library(corrplot)
library(mice)
library(pastecs)
library(knitr)
library(mvtnorm)
library(MASS)

library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

library(SimDesign)
library(gfoRmula)

rm(list = ls())

generate_panel_data <- function(
    popN = 1e6,
    T = 5,
    outcome_type = "binary",               # "binary" or "continuous"
    relationship_type = "quadratic",      # "linear", "quadratic", "quadratic+interaction"
    confounding = "high"                  # "low" or "high"
) {
  library(MASS)
  library(dplyr)
  library(corrplot)
  
  p <- 3
  k <- popN
  
  # ---- Choose confounding level via exposure correlation ----------------
  covmat <- switch(confounding,
                   "low"  = toeplitz(c(1, .3, .2)),
                   "high" = matrix(c(1,.05, .2, .05, 1, .1, .2, .1,1), nrow = 3, byrow= T))
  
  sigsq.trueL <- 1
  sigsq.trueY <- 1
  sigsq.trueA <- rep(0.05, T - 1)
  
  set.seed(1)
  dat <- list()
  dat$sex <- rbinom(k, 1, 0.5)
  dat$C0 <- rnorm(k, mean = 0.25 * (1 - dat$sex), sd = 1)
  dat$At0 <- mvtnorm::rmvnorm(k, rep(0, p), sigma = covmat)
  colnames(dat$At0) <- paste0("At0_", 1:p)
  
  make_hfunL <- function(type) {
    if (type == "linear") {
      function(z) 0.25 * z[1] + 0.5 * z[2]
    } else if (type == "quadratic") {
      function(z) 0.25 * z[1] + 0.5 * z[2] + 0.25 * z[2]^2
    } else if (type == "quadratic+interaction") {
      function(z) 0.25 * z[1] + 0.5 * z[2] + 0.25 * z[1] * z[2] + 0.25 * z[2]^2
    }
  }
  
  for (t in 1:(T - 1)) {
    prev_At <- dat[[paste0("At", t - 1)]]
    prev_L  <- if (t == 1) dat$C0 else dat[[paste0("L", t - 1)]]
    hAt_t <- matrix(NA, nrow = k, ncol = p)
    for (j in 1:p) {
      hAt_t[, j] <- 1/3 * prev_At[, j]^2 + 0.1 * prev_At[, j] +
        0.1 * dat$C0 + 0.1 * prev_L
    }
    dat[[paste0("hAt", t)]] <- hAt_t
    dat[[paste0("epsAt", t)]] <- matrix(rnorm(k * p, sd = sqrt(sigsq.trueA[t])), k, p)
    dat[[paste0("At", t)]] <- hAt_t + dat[[paste0("epsAt", t)]]
    colnames(dat[[paste0("At", t)]]) <- paste0("At", t, "_", 1:p)
    
    hfunL <- make_hfunL(relationship_type)
    dat[[paste0("hL", t)]] <- apply(cbind(dat[[paste0("At", t)]][, 2], prev_L), 1, hfunL)
    dat[[paste0("L", t)]] <- dat[[paste0("hL", t)]] + rnorm(k, sd = sqrt(sigsq.trueL))
  }
  
  At_list <- lapply(0:(T - 1), function(t) dat[[paste0("At", t)]])
  names(At_list) <- paste0("At", 0:(T - 1))
  L_list <- c(list(L0 = dat$C0), dat[paste0("L", 1:(T - 1))])
  names(L_list) <- paste0("L", 0:(T - 1))
  
  dat$ALL <- do.call(cbind, c(At_list, L_list, list(dat$sex)))
  colnames(dat$ALL) <- c(
    unlist(lapply(0:(T - 1), function(t) paste0("At", t, "_", 1:p))),
    paste0("L", 0:(T - 1)),
    "sex"
  )
  
  hfunY <- function(z) {
    out <- 0
    for (t in 0:(T - 1)) {
      L <- z[[paste0("L", t)]]
      M2 <- z[[paste0("At", t, "_2")]]
      if (relationship_type == "linear") {
        out <- out + 0.25 * M2 + 0.25 * L
      } else if (relationship_type == "quadratic") {
        out <- out + 0.25 * M2 + 0.25 * L + 0.25 * L^2
      } else if (relationship_type == "quadratic+interaction") {
        out <- out + 0.25 * M2 + 0.25 * L + 0.25 * L^2 + 0.1 * M2 * L
      }
    }
    return(out)
  }
  
  dat$hY <- apply(dat$ALL, 1, function(row) hfunY(as.list(row)))
  dat$epsY <- rnorm(k, sd = sqrt(sigsq.trueY))
  dat$z <- dat$hY + dat$epsY
  
  dat$y <- if (outcome_type == "binary") {
    rbinom(k, 1, pnorm(dat$z))
  } else {
    dat$z
  }
  
  for (t in 0:(T - 1)) {
    colnames(dat[[paste0("At", t)]]) <- paste0("logM", 1:p, "_", t)
  }
  names(L_list) <- c("waist0", paste0("waist", 1:(T - 1)))
  
  df <- data.frame(
    sex = dat$sex,
    waist0 = dat$C0,
    do.call(cbind, dat[paste0("At", 0:(T - 1))]),
    do.call(cbind, L_list[-1]),
    Y = dat$y
  )
  df$id <- 1:k
  
  colnames(df) <- c("sex", "waist0",
                    unlist(lapply(0:(T - 1), function(t) paste0("logM", 1:p, "_", t))),
                    paste0("waist", 1:(T - 1)), "Y", "id")
  
  corrplot(cor(df[, -ncol(df)]))
  saveRDS(df, file = paste0("popn_", T, "t_", relationship_type, "_", confounding, "_", outcome_type, ".rds"))
  
  return(df)
}

# Generate a binary outcome dataset with low confounding, 5 time points, quadratic+interaction relationship
df <- generate_panel_data(
  popN = 1e6,
  T = 3,
  outcome_type = "continuous",
  relationship_type = "quadratic+interaction",
  confounding = "high"
)

