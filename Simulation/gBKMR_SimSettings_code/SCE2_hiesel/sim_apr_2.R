args<-commandArgs(TRUE)
currind <-as.integer(args[1])
print(currind)

library(bkmr)
library(data.table)
library(CBPS)
library(dplyr)

sim_popn = readRDS("popn_lowcor_nonlinearLY_strongL_1.rds")
n = 500

sel<-seq(10000,12000,by=25)

sim_res = rep(NA, 3)
  
set.seed(currind)

dat_sim = sim_popn[sample(sim_popn$id, n, replace=F),]  
 
 
Y = dat_sim$Y
Cov_mat_y = dplyr::select(dat_sim,sex,  waist0 )
#Exp_mat_y_withL= dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, 
#                         logM1_1, logM2_1, logM3_1 , waist1  )
Exp_mat_y_noL= dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, 
                               logM1_1, logM2_1, logM3_1 )



#fitkm_y_withL <- kmbayes(y = Y, Z = Exp_mat_y_withL, X = Cov_mat_y, iter = 12000, verbose = FALSE, 
#                     groups = c(1, 1, 1, 2, 2, 2, 2),varsel = TRUE) 


fitkm_y_noL <- kmbayes(y = Y, Z = Exp_mat_y_noL, X = Cov_mat_y, iter = 12000, verbose = FALSE, 
                         groups = c(1, 1, 1, 2, 2, 2),varsel = TRUE) 

#saveRDS(fitkm_y_withL, file = paste0("fitkm_y_withL_", currind, ".rds")) 
saveRDS(fitkm_y_noL, file = paste0("fitkm_y_noL_", currind, ".rds")) 
 
#fitkm_y <- readRDS(paste0("fitkm_y_", currind, ".rds"))

 start.time <- proc.time()
K = 1000

 A <- dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0,  logM1_1, logM2_1, logM3_1)

a <-   apply(A, 2, quantile, probs=0.25)
astar <-   apply(A, 2, quantile, probs=0.75)
 
X.predict.Y <- matrix(colMeans(Cov_mat_y),nrow=1)

newz      <- rbind(a,astar)
set.seed(93020)
Y.samp <- SamplePred(fitkm_y_noL, Znew =newz, Xnew = X.predict.Y, sel=sel)
Ya       <- as.vector(Y.samp[,"znew1"])
Yastar   <- as.vector(Y.samp[,"znew2"])

 
saveRDS(Ya, paste0("Ya_", currind, ".rds"))
saveRDS(Yastar, paste0("Yastar_", currind, ".rds"))






