args<-commandArgs(TRUE)
currind <-as.integer(args[1])
print(currind)

library(bkmr)
library(data.table)
library(CBPS)
library(dplyr)

sim_popn = readRDS("popn_highcor_nonlinearLY_3.rds")
n = 500

sel<-seq(10000,12000,by=25)

sim_res = rep(NA, 3)


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

para_l1 = apply(fitkm_l1$beta[6000:12000,],2,mean)
#para_l2 = apply(fitkm_l2$beta[40000:50000,],2,mean)
para_y = apply(fitkm_y$beta[6000:12000,],2,mean)
para_l1l2y = c( para_l1,  para_y)

write.csv(para_l1l2y, file = paste0("para_res", currind, ".csv"))


#####g-formula

library(gfoRmula)
simdat_long <- reshape(dat_sim, direction='long',
                       varying=list(waist = c('waist0', 'waist1' ),
                                    logM1 = c('logM1_0', 'logM1_1'  ),
                                    logM2 = c('logM2_0', 'logM2_1'  ),
                                    logM3 = c('logM3_0', 'logM3_1'  )),
                       timevar='time',
                       times=c('0', '1' ),
                       v.names=c('waist', 'logM1','logM2' ,'logM3'  ),
                       idvar='id')

simdat_long$time = as.numeric(simdat_long$time)
simdat_long= simdat_long[order(simdat_long[,3], simdat_long[,4]), ]
simdat_long = as.data.table(simdat_long)
#glimpse(simdat_long)

id <- 'id'
time_name <- 'time'
covnames <- c('waist','logM1', 'logM2', 'logM3')
outcome_name <- 'Y'

histories <- c(lagged)
histvars <- list(c('waist','logM1', 'logM2', 'logM3'))

covtypes <- c('normal', 'normal', 'normal', 'normal')
covparams <- list(covmodels = c(waist ~ lag1_logM1 + lag1_logM2 +lag1_logM3 +  lag1_waist +
                                  sex + time,
                                logM1 ~ lag1_logM1 + lag1_logM2 +lag1_logM3 +
                                  waist+ lag1_waist+ sex + time,
                                logM2 ~ lag1_logM1 + lag1_logM2 +lag1_logM3 +
                                  waist+ lag1_waist+ sex+  time,
                                logM3 ~ lag1_logM1 + lag1_logM2 +lag1_logM3 +
                                  waist+ lag1_waist+ sex+  time))

ymodel <- Y ~ logM1 + logM2 + logM3 + waist + lag1_waist +lag1_logM1 + lag1_logM2 +lag1_logM3+ sex



Ms = dplyr::select(simdat_long, logM1, logM2, logM3)
quants = apply(Ms , 2 , quantile , probs = c(.25, .75) , na.rm = TRUE )
M25 <-  quants[1, ]
M75 <-  quants[2, ]

intvars <- list(c('logM1', 'logM2', 'logM3'), c('logM1', 'logM2', 'logM3'))
interventions <- list(list(c(threshold, M25, M25)), list(c(threshold, M75, M75)))
int_descript <- c('threshold-first quantile', 'threshold-third quantile')

nsimul <- 10000
gform_cont_eof <- gformula_continuous_eof(obs_data = simdat_long, id = id,
                                          time_name = time_name, covnames = covnames,
                                          outcome_name = outcome_name, covtypes = covtypes,
                                          covparams = covparams, ymodel = ymodel,
                                          intvars = intvars, interventions = interventions,
                                          int_descript = int_descript, histories = histories,
                                          histvars = histvars, basecovs = c("sex"),
                                          seed = 1234, parallel = FALSE, nsamples = 5,
                                          nsimul = nsimul )


# gform_cont_eof$result

g_yhat_25 = as.numeric( gform_cont_eof$result[2,4])
g_yhat_75 = as.numeric( gform_cont_eof$result[3,4])


# gform_cont_eof$result

g_yhat_25 = as.numeric( gform_cont_eof$result[2,4])
g_yhat_75 = as.numeric( gform_cont_eof$result[3,4])
diff_gform = g_yhat_75 - g_yhat_25


##### Linear Regression


# df_lm = dplyr::mutate(dat_sim,
#                       glucose0 = glucose0 - mean(glucose0),
#                        waist0 =  waist0 - mean(waist0), waist1 = waist1 - mean(waist1) ,
#                        age = age - mean(age))



lm_model = lm(Y ~  as.factor(sex)   +waist0 +waist1 +
                logM1_0+logM2_0+logM3_0+logM1_1+logM2_1+logM3_1 , dat = dat_sim)
summary(lm_model)
newdat_25 <- data.frame(  sex = as.factor(0),
                        logM1_0 = quantile(dat_sim$logM1_0)[2],
                        logM1_1 = quantile(dat_sim$logM1_1)[2],
                        logM2_0 = quantile(dat_sim$logM2_0)[2],
                        logM2_1 = quantile(dat_sim$logM2_1)[2],
                        logM3_0 = quantile(dat_sim$logM3_0)[2],
                        logM3_1 = quantile(dat_sim$logM3_1)[2],
                        waist0 = mean(dat_sim$waist0),
                        waist1 = mean(dat_sim$waist1)
)

newdat_75 <- data.frame( sex = as.factor(0),
                        logM1_0 = quantile(dat_sim$logM1_0)[4],
                        logM1_1 = quantile(dat_sim$logM1_1)[4],
                        logM2_0 = quantile(dat_sim$logM2_0)[4],
                        logM2_1 = quantile(dat_sim$logM2_1)[4],
                        logM3_0 = quantile(dat_sim$logM3_0)[4],
                        logM3_1 = quantile(dat_sim$logM3_1)[4],
                        waist0 = mean(dat_sim$waist0),
                        waist1 = mean(dat_sim$waist1) )
yhat_lm_25 = predict(lm_model, newdat_25 )
yhat_lm_75 = predict(lm_model, newdat_75 )
diff_lm = yhat_lm_75-yhat_lm_25





sim_res[1] <-  diff_gBKMR
sim_res[2] <- diff_gform
sim_res[3] <- diff_lm

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









