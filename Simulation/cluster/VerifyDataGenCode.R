library(bkmr)
library(data.table)
library(CBPS)
library(dplyr)
library(rstudioapi)
library(causalbkmr)

rm(list = ls())

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

nSim <- 100

############################################
# 1) VERIFY SIMULATED DATA
############################################
popn <- readRDS("popn_3t_highconf_linearLY_binary_probit.rds")

summary(popn$waist0)  # Should be ~ N(90, 15^2)
table(popn$sex)       # Expect about half 0, half 1
#check correlations among baseline exposures (roughly no correlation if you used diag(3)):
cor(popn[, c("logM1_0", "logM2_0", "logM3_0")])
#quick plot of a relevant relationship
plot(popn$waist0, popn$waist1,
     main = "waist1 vs waist0",
     xlab = "waist0", ylab = "waist1")

plot(popn$waist1, popn$waist2,
     main = "waist2 vs waist1",
     xlab = "waist1", ylab = "waist2")

############################################
# 2) LOOP OVER THE 100 REPLICATES AND LOAD RESULTS
############################################

# Create lists (or vectors) to store results
fitkm_l1_list      <- vector("list", nSim)
fitkm_l2_list      <- vector("list", nSim)
fitkm_y_list       <- vector("list", nSim)
para_res_list      <- vector("list", nSim)
risks_l1_list      <- vector("list", nSim)
risks_y_list       <- vector("list", nSim)
sim_res_list       <- vector("list", nSim)
YaLa_list          <- vector("list", nSim)
YastarLastar_list  <- vector("list", nSim)

for(i in seq_len(nSim)) {
  # Adjust filenames if needed; this is just a template
  fitkm_l1_list[[i]]     <- readRDS(paste0("s6_result/fitkm_l1_", i, ".rds"))
  fitkm_l2_list[[i]]     <- readRDS(paste0("s6_result/fitkm_l2_", i, ".rds"))
  fitkm_y_list[[i]]      <- readRDS(paste0("s6_result/fitkm_y_", i, ".rds"))
  para_res_list[[i]]     <- read.csv(paste0("s6_result/para_res", i, ".csv"))
  risks_l1_list[[i]]     <- readRDS(paste0("s6_result/risks_singvar_l1_", i, ".rds"))
  risks_y_list[[i]]      <- readRDS(paste0("s6_result/risks_singvar_y_", i, ".rds"))
  sim_res_list[[i]]      <- read.csv(paste0("s6_result/Sim_res", i, ".csv"))
  YaLa_list[[i]]         <- readRDS(paste0("s6_result/YaLa_", i, ".rds"))
  YastarLastar_list[[i]] <- readRDS(paste0("s6_result/YastarLastar_", i, ".rds"))
}

############################################
# 3) CHECK MODEL FITS & INFERENCES
############################################

# 3a) Example: Summaries of the BKMR fits in a single replicate
#    (Do this for a few replicates or in a loop if you prefer.)
summary(fitkm_l1_list[[1]])
summary(fitkm_l2_list[[1]])
summary(fitkm_y_list[[1]])

# 3b) Check parameter estimates across simulations

all_para_res <- dplyr::bind_rows(para_res_list, .id = "sim_id")
summary(all_para_res)

all_sim_res <- dplyr::bind_rows(sim_res_list, .id = "sim_id")
summary(all_sim_res)

# If these CSVs contain metrics like bias, coverage, MSE, etc., you can summarize:
# e.g., mean bias or coverage across the 100 replications:
all_sim_res %>%
  group_by(parameter_name) %>%
  summarize(mean_est = mean(estimate, na.rm = TRUE),
            sd_est   = sd(estimate, na.rm = TRUE),
            mean_bias = mean(bias, na.rm = TRUE),
            coverage  = mean(coverage, na.rm = TRUE))

# 3c) Examine single-variable risk summaries from the BKMR fits:
# e.g., "risks_singvar_l1_" might have predicted risk changes for L1
# You can check if the sign, magnitude, or rank ordering is aligned with the known data-generating process
# (depending on what "risks_singvar_l1_" actually stores).



