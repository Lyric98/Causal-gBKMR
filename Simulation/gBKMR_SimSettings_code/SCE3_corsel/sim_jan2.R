args<-commandArgs(TRUE)
currind <-as.integer(args[1])
print(currind)

library(bkmr)
library(data.table)
library(CBPS)
library(dplyr)

sim_popn = readRDS("popn_3t_highconf_quadraticLY.rds")
n = 500

sel<-seq(22000,24000,by=25)

sim_res = rep(NA, 2)

set.seed(currind)
dat_sim = sim_popn[sample(sim_popn$id, n, replace=F),]


L1 = dat_sim$waist1

Cov_mat_l1 = dplyr::select(dat_sim,   sex,  waist0 )
Exp_mat_l1 = dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0 )



L2 = dat_sim$waist1

Cov_mat_l2 = dplyr::select(dat_sim,   sex,  waist0 )
Exp_mat_l2 = dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, logM1_1, logM2_1, logM3_1, waist1 )


Y = dat_sim$Y
Cov_mat_y = dplyr::select(dat_sim,  sex,  waist0 )
Exp_mat_y= dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, logM1_1, logM2_1, logM3_1 ,
                         logM1_2, logM2_2, logM3_2 , waist1, waist2  )


  fitkm_l1 <- kmbayes(y = L1, Z = Exp_mat_l1, X = Cov_mat_l1, iter = 24000, verbose = FALSE, varsel = TRUE)
 saveRDS(fitkm_l1, file = paste0("fitkm_l1_", currind, ".rds"))
 #
 #
 fitkm_l2 <- kmbayes(y = L2, Z = Exp_mat_l2, X = Cov_mat_l2, iter = 24000, verbose = FALSE, varsel = TRUE)
 saveRDS(fitkm_l2, file = paste0("fitkm_l2_", currind, ".rds"))

  fitkm_y <- kmbayes(y = Y, Z = Exp_mat_y, X = Cov_mat_y, iter = 24000, verbose = FALSE, varsel = TRUE)
 saveRDS(fitkm_y, file = paste0("fitkm_y_", currind, ".rds"))

#fitkm_l1 <- readRDS(paste0("fitkm_l1_", currind, ".rds"))
# fitkm_l2 <- readRDS(paste0("fitkm_l2_", currind, ".rds"))
#fitkm_y <- readRDS(paste0("fitkm_y_", currind, ".rds"))

  start.time <- proc.time()
K = 1000

 A <- dplyr::select(dat_sim, logM1_0, logM2_0, logM3_0, logM1_1, logM2_1, logM3_1, logM1_2, logM2_2, logM3_2)

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


L2aL1a.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)
L2astarL1astar.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)


for(j in 1:length(sel)){
  print(j)
  L1a.j     <- L1a.samp[j,]
  L1astar.j <- L1astar.samp[j,]

  aL1a.j         <- cbind(matrix(a[1:6], nrow=K, ncol=length(a[1:6]), byrow=TRUE), L1a.j)
  astarL1astar.j <- cbind(matrix(astar[1:6], nrow=K, ncol=length(astar[1:6]), byrow=TRUE), L1astar.j)

  for(k in 1:K){
    newz <- rbind(aL1a.j[k,], astarL1astar.j[k,])

    set.seed(j + 10000)
    L2astarL1astar.jk <- SamplePred(fitkm_l2, Znew = newz, Xnew = X.predict.L, sel=sel[j])

    L2aL1a.samp.mat[j,k]<- L2astarL1astar.jk[,"znew1"]
    L2astarL1astar.samp.mat[j,k]<- L2astarL1astar.jk[,"znew2"]

  }
  end.time.temp <- proc.time()
  if(j%%50==0) print(paste("iter", j, "time: ", round((end.time.temp - start.time)["elapsed"]/60,2),"min"))
}


sigma.samp.l2  <- sqrt(fitkm_l2$sigsq.eps[sel])

L2aL1a = as.vector(apply(L2aL1a.samp.mat, 1, mean))
L2astarL1astar = as.vector(apply(L2astarL1astar.samp.mat, 1, mean))

random.samp.l2.1 <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)
random.samp.l2.2 <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)

L2aL1a.samp  <- L2aL1a+ sigma.samp.l2*random.samp.l2.1
L2astarL1astar.samp  <- L2astarL1astar+ sigma.samp.l2*random.samp.l2.2


YaLa.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)
YastarLastar.samp.mat <- matrix(NA,nrow=length(sel),ncol=K)


for(j in 1:length(sel)){
  print(j)
  L1a.j     <- L1a.samp[j,]
  L2a.j     <- L2aL1a.samp[j,]
  aLa.j     <- cbind(matrix(a, nrow=K, ncol=length(a), byrow=TRUE), L1a.j, L2a.j)


  L1astar.j     <- L1astar.samp[j,]
  L2astar.j     <- L2astarL1astar.samp[j,]
  astarLastar.j     <- cbind(matrix(astar, nrow=K, ncol=length(astar), byrow=TRUE), L1astar.j, L2astar.j)


  for(k in 1:K){
    newz <- rbind(aLa.j[k,], astarLastar.j[k,])

    set.seed(j + 10000)
    Y.jk <- SamplePred(fitkm_y, Znew = newz, Xnew = X.predict.Y, sel=sel[j])

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

para_l1 = apply(fitkm_l1$beta[22000:24000,],2,mean)
para_l2 = apply(fitkm_l2$beta[22000:24000,],2,mean)
para_y = apply(fitkm_y$beta[22000:24000,],2,mean)
para_l1l2y = c( para_l1,  para_l2, para_y)

write.csv(para_l1l2y, file = paste0("para_res", currind, ".csv"))


##### Linear Regression


# df_lm = dplyr::mutate(dat_sim,
#                       glucose0 = glucose0 - mean(glucose0),
#                        waist0 =  waist0 - mean(waist0), waist1 = waist1 - mean(waist1) ,
#                        age = age - mean(age))



lm_model = lm(Y ~  as.factor(sex)   +waist0 +waist1 +waist2+
                logM1_0+logM2_0+logM3_0+logM1_1+logM2_1+logM3_1+logM1_2+logM2_2+logM3_2 , dat = dat_sim)
#summary(lm_model)
newdat_25 <- data.frame(  sex = as.factor(0),
                          logM1_0 = quantile(dat_sim$logM1_0)[2],
                          logM1_1 = quantile(dat_sim$logM1_1)[2],
                          logM1_2 = quantile(dat_sim$logM1_2)[2],
                          logM2_0 = quantile(dat_sim$logM2_0)[2],
                          logM2_1 = quantile(dat_sim$logM2_1)[2],
                          logM2_2 = quantile(dat_sim$logM2_2)[2],
                          logM3_0 = quantile(dat_sim$logM3_0)[2],
                          logM3_1 = quantile(dat_sim$logM3_1)[2],
                          logM3_2 = quantile(dat_sim$logM3_2)[2],
                        waist0 = mean(dat_sim$waist0),
                        waist1 = mean(dat_sim$waist1),
                        waist2 = mean(dat_sim$waist2)
)

newdat_75 <- data.frame( sex = as.factor(0),
                         logM1_0 = quantile(dat_sim$logM1_0)[4],
                         logM1_1 = quantile(dat_sim$logM1_1)[4],
                         logM1_2 = quantile(dat_sim$logM1_1)[4],
                         logM2_0 = quantile(dat_sim$logM2_0)[4],
                         logM2_1 = quantile(dat_sim$logM2_1)[4],
                         logM2_2 = quantile(dat_sim$logM2_2)[4],
                         logM3_0 = quantile(dat_sim$logM3_0)[4],
                         logM3_1 = quantile(dat_sim$logM3_1)[4],
                         logM3_2 = quantile(dat_sim$logM3_2)[4],
                        waist0 = mean(dat_sim$waist0),
                        waist1 = mean(dat_sim$waist1) ,
                        waist2 = mean(dat_sim$waist2)
                        )
yhat_lm_25 = predict(lm_model, newdat_25 )
yhat_lm_75 = predict(lm_model, newdat_75 )
diff_lm = yhat_lm_75-yhat_lm_25





sim_res[1] <-  diff_gBKMR
sim_res[2] <- diff_lm

write.csv(sim_res, file = paste0("Sim_res", currind, ".csv"))




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









