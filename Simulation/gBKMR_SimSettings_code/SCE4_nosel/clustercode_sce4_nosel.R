args<-commandArgs(TRUE)
currind <-as.integer(args[1])
print(currind)

library(bkmr)
library(data.table)
library(CBPS)
library(dplyr)

sim_popn = readRDS("popn_highcor_nonlinearLY_3_binary_probit.rds")
n = 500

sel<-seq(10000,12000,by=25)
 
 
set.seed(currind)
dat_sim = sim_popn[sample(sim_popn$id, n, replace=F),]  


L1 = dat_sim$waist1

Cov_mat_l1 = dplyr::select(dat_sim,   sex,  waist0 )
Exp_mat_l1 = dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0 )

 

Y = dat_sim$Y
Cov_mat_y = dplyr::select(dat_sim,sex,  waist0 )
Exp_mat_y= dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, 
                         logM1_1, logM2_1, logM3_1 , waist1  )


  fitkm_l1 <- kmbayes(y = L1, Z = Exp_mat_l1, X = Cov_mat_l1, iter = 12000, verbose = FALSE, varsel = FALSE) 
 saveRDS(fitkm_l1, file = paste0("fitkm_l1_", currind, ".rds"))
 #  
 # 
  fitkm_y <- kmbayes(y = Y, Z = Exp_mat_y, X = Cov_mat_y, iter = 12000, verbose = FALSE, varsel = FALSE) 
 saveRDS(fitkm_y, file = paste0("fitkm_y_", currind, ".rds")) 

#fitkm_l1 <- readRDS(paste0("fitkm_l1_", currind, ".rds"))
# fitkm_l2 <- readRDS(paste0("fitkm_l2_", currind, ".rds"))
#fitkm_y <- readRDS(paste0("fitkm_y_", currind, ".rds"))

  start.time <- proc.time()
K = 1000

 A <- dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0,  logM1_1, logM2_1, logM3_1)

a <-   apply(A, 2, quantile, probs=0.25)
astar <-   apply(A, 2, quantile, probs=0.75)

# Calculate L_a*

X.predict.L <- matrix(colMeans(Cov_mat_l1),nrow=1)
X.predict.Y <- matrix(colMeans(Cov_mat_y),nrow=1)

newz      <- rbind(a,astar)
set.seed(93020)
EL1.samp <- SamplePred(fitkm_l1, Znew =newz[,1:3], Xnew = X.predict.L, sel=sel)
L1a       <- as.vector(EL1.samp[,"znew1"])
L1astar   <- as.vector(EL1.samp[,"znew2"])

sigma.samp   <- sqrt(fitkm_l1$sigsq.eps[sel])
set.seed(93020)
random.samp1 <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)
random.samp2 <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)

L1a.samp     <- L1a + sigma.samp*random.samp1
L1astar.samp <- L1astar + sigma.samp*random.samp2

YaLa.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)
YastarLastar.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)
  

for(j in 1:length(sel)){
  print(j)
  L1a.j     <- L1a.samp[j,]
  L1astar.j <- L1astar.samp[j,]
  
  aL1a.j         <- cbind(matrix(a, nrow=K, ncol=length(a), byrow=TRUE), L1a.j) 
  astarL1astar.j <- cbind(matrix(astar, nrow=K, ncol=length(astar), byrow=TRUE), L1astar.j) 
  
  for(k in 1:K){
    newz <- rbind(aL1a.j[k,], astarL1astar.j[k,])
    
    set.seed(j + 10000)
    Y.jk <- SamplePred(fitkm_y, Znew = newz, Xnew = X.predict.Y, sel=sel[j], type = "response")
    
    YaLa.samp.mat[j,k]<- Y.jk[,"znew1"]
    YastarLastar.samp.mat[j,k]<- Y.jk[,"znew2"] 
    
  }   
  end.time.temp <- proc.time() 
  if(j%%50==0) print(paste("iter", j, "time: ", round((end.time.temp - start.time)["elapsed"]/60,2),"min")) 
}


YaLa         <- apply(YaLa.samp.mat,        1,mean)
YastarLastar <- apply(YastarLastar.samp.mat,1,mean)

saveRDS( YaLa, paste0("YaLa_", currind, ".rds"))
saveRDS( YastarLastar, paste0("YastarLastar_", currind, ".rds"))

##g-BKMR
diff_gBKMR = mean(YastarLastar) - mean(YaLa) 

para_l1 = apply(fitkm_l1$beta[6000:12000,],2,mean)
#para_l2 = apply(fitkm_l2$beta[40000:50000,],2,mean)
para_y = apply(fitkm_y$beta[6000:12000,],2,mean)
para_l1l2y = c( para_l1,  para_y)

write.csv(para_l1l2y, file = paste0("para_res", currind, ".csv")) 
 

risks_singvar_l1 <- SingVarRiskSummaries(fit = fitkm_l1, y = L1, Z = Exp_mat_l1, X = as.matrix(Cov_mat_l1), 
                                         qs.diff = c(0.25,0.75), 
                                         q.fixed = c( 0.25, 0.5, 0.75),
                                         method = "exact")
saveRDS(risks_singvar_l1, file = paste0("risks_singvar_l1_", currind, ".rds")) 


risks_singvar_y <- SingVarRiskSummaries(fit = fitkm_y, y = Y, Z = Exp_mat_y, X = as.matrix(Cov_mat_y), 
                                        qs.diff = c(0.25,0.75), 
                                        q.fixed = c( 0.25,0.5, 0.75),
                                        method = "exact")

saveRDS(risks_singvar_y, file = paste0("risks_singvar_y_", currind, ".rds")) 









