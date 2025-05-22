run_gbkmr_panel <- function(
    sim_popn,
    T = 5,
    currind = 1,
    sel = seq(22000, 24000, by = 25),
    n = 500,
    K = 1000,
    iter = 24000,
    parallel = TRUE,
    save_exposure_preds = TRUE,
    return_ci = TRUE,
    make_plots = TRUE,
    use_knots = TRUE,
    n_knots = 50
) {
  library(bkmr)
  library(fields)
  library(dplyr)
  library(ggplot2)
  library(parallel)
  
  if (max(sel) > iter) stop("sel contains indices beyond total MCMC iterations!")
  
  message("Subsampling population...")
  set.seed(currind)
  dat_sim <- sim_popn[sample(sim_popn$id, n, replace = FALSE), ]
  cov_names_common <- c("sex", "waist0")
  X_common <- dplyr::select(dat_sim, all_of(cov_names_common))
  
  fitkm_list <- list()
  L_values_a <- list()
  L_values_astar <- list()
  
  for (t in 1:(T - 1)) {
    message(paste("Fitting mediator model L", t))
    y_L <- dat_sim[[paste0("waist", t)]]
    exp_names <- unlist(lapply(0:(t - 1), function(s) paste0("logM", 1:3, "_", s)))
    if (t > 1) exp_names <- c(exp_names, paste0("waist", 1:(t - 1)))
    Z <- dplyr::select(dat_sim, all_of(exp_names))
    X <- X_common
    
    scale_Z <- scale(Z)
    attr_list <- list(center = attr(scale_Z, "scaled:center"),
                      scale = attr(scale_Z, "scaled:scale"))
    knots <- if (use_knots) {
      set.seed(1000)
      fields::cover.design(scale_Z, nd = n_knots)$design
    } else NULL
    
    fit <- kmbayes(y = y_L, Z = scale_Z, X = X, iter = iter, varsel = TRUE,
                   verbose = FALSE, knots = knots)
    fitkm_list[[paste0("L", t)]] <- list(fit = fit, scale_info = attr_list)
  }
  
  message("Fitting outcome model Y")
  Y <- dat_sim$Y
  exp_names_y <- c(unlist(lapply(0:(T - 1), function(s) paste0("logM", 1:3, "_", s))),
                   paste0("waist", 1:(T - 1)))
  Z_y <- dplyr::select(dat_sim, all_of(exp_names_y))
  scale_Zy <- scale(Z_y)
  scale_info_y <- list(center = attr(scale_Zy, "scaled:center"),
                       scale = attr(scale_Zy, "scaled:scale"))
  knots_y <- if (use_knots) {
    set.seed(1000)
    fields::cover.design(scale_Zy, nd = n_knots * 2)$design
  } else NULL
  
  fit_y <- kmbayes(y = Y, Z = scale_Zy, X = X_common, iter = iter, varsel = TRUE,
                   verbose = FALSE, knots = knots_y)
  fitkm_y <- list(fit = fit_y, scale_info = scale_info_y)
  
  # Counterfactuals
  message("Computing counterfactual exposure vectors...")
  A <- dplyr::select(dat_sim, starts_with("logM"))
  a <- apply(A, 2, quantile, probs = 0.25)
  astar <- apply(A, 2, quantile, probs = 0.75)
  
  X.predict <- matrix(colMeans(X_common), nrow = 1)
  
  for (t in 1:(T - 1)) {
    message(paste("Predicting mediator L", t))
    scale_info <- fitkm_list[[paste0("L", t)]]$scale_info
    fit <- fitkm_list[[paste0("L", t)]]$fit
    
    if (t == 1) {
      newz <- rbind(a[1:3], astar[1:3])
    } else {
      prev_La <- sapply(1:(t - 1), function(j) L_values_a[[j]])
      prev_Lastar <- sapply(1:(t - 1), function(j) L_values_astar[[j]])
      newz <- rbind(
        c(a[1:(3 * t)], rowMeans(t(prev_La))),
        c(astar[1:(3 * t)], rowMeans(t(prev_Lastar)))
      )
    }
    
    Znew_scaled <- scale(newz, center = scale_info$center, scale = scale_info$scale)
    L_pred <- SamplePred(fit, Znew = Znew_scaled, Xnew = X.predict, sel = sel)
    L_values_a[[t]] <- as.vector(L_pred[, "znew1"])
    L_values_astar[[t]] <- as.vector(L_pred[, "znew2"])
  }
  
  message("Predicting outcome Y")
  Z_final <- rbind(
    c(a, sapply(L_values_a, mean)),
    c(astar, sapply(L_values_astar, mean))
  )
  Z_final_scaled <- scale(Z_final, center = fitkm_y$scale_info$center,
                          scale = fitkm_y$scale_info$scale)
  Y_pred <- SamplePred(fitkm_y$fit, Znew = Z_final_scaled, Xnew = X.predict, sel = sel)
  
  Ya <- Y_pred[, "znew1"]
  Yastar <- Y_pred[, "znew2"]
  diff_gBKMR <- mean(Yastar) - mean(Ya)
  
  message(sprintf("g-BKMR effect estimate: %.4f", diff_gBKMR))
  
  if (make_plots) {
    df_plot <- data.frame(Scenario = c("a", "astar"), Mean = c(mean(Ya), mean(Yastar)))
    print(ggplot(df_plot, aes(x = Scenario, y = Mean)) +
            geom_col(fill = "skyblue") +
            theme_minimal() +
            ggtitle("Counterfactual Means"))
  }
  
  return(list(
    diff_gBKMR = diff_gBKMR,
    Ya = Ya,
    Yastar = Yastar,
    beta_all = c(
      unlist(lapply(fitkm_list, function(l) colMeans(l$fit$beta[sel, , drop = FALSE]))),
      colMeans(fitkm_y$fit$beta[sel, , drop = FALSE])
    ),
    L_values_a = L_values_a,
    L_values_astar = L_values_astar
  ))
}
