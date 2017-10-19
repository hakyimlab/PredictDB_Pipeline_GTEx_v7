p_vals <- numeric()
for (i in 1:5000) {
  x <- rnorm(100)
  y <- rnorm(100)
  p_val <- summary(lm(y ~ x))$coef[2,4]
  p_vals <- c(p_vals, p_val)
}

hist(p_vals, breaks = 20, main= 'Ordinary Least Squares - Null')

library(glmnet)

evaluate_performance <- function(cis_gt, expr_adj, fit, best_lam_ind, cv_fold_ids, n_folds) {
  n_nonzero <- fit$nzero[best_lam_ind]
  if (n_nonzero > 0) {
    R2 <- rep(0, n_folds)
    for (j in (1:n_folds)) {
      fold_idxs <- which(cv_fold_ids[,1] == j)
      tss <- sum(expr_adj[fold_idxs]**2)
      rss <- sum((expr_adj[fold_idxs] - fit$fit.preval[fold_idxs, best_lam_ind])**2)
      R2[j] <- 1 - (rss/tss)
    }
    best_fit <- fit$glmnet.fit
    expr_adj_pred <- predict(best_fit, cis_gt, s = best_lambda)
    tss <- sum((expr_adj-mean(expr_adj))**2)
    rss <- sum((expr_adj - expr_adj_pred)**2)
    
    n_samp <- length(expr_adj)
    #f_stat <- ((tss - rss) / n_nonzero) / (rss / (n_samp - n_nonzero - 1))
    #p_val <- pf(f_stat, n_samp - 1, n_samp - n_nonzero - 1, lower.tail = FALSE)
    weights <- best_fit$beta[which(best_fit$best[,best_lam_ind] != 0), best_lam_ind]
    R2_mean <- mean(R2)
    R2_sd <- sd(R2)
    inR2 <- 1 - (rss/tss)
    # Old way
    pred_perf <- summary(lm(expr_adj ~ fit$fit.preval[,best_lam_ind]))
    pred_perf_rsq <- pred_perf$r.squared
    
    pred_perf_pval <- pred_perf$coef[2,4]
    #one_sided_pval <- cor.test(expr_adj, fit$fit.preval[,best_lam_ind], alternative = 'greater')$p.value
    out <- list(weights = weights, n_weights = n_nonzero, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pred_perf_rsq = pred_perf_rsq, pred_perf_pval = pred_perf_pval)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pred_perf_rsq = NA, pred_perf_pval = NA)
  }
  out
}

n_samples <- 150
n_times <- 1
cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
n_k_folds <- 10
for (j in 1:n_times)
  cv_fold_ids[,j] <- sample(1:n_k_folds, n_samples, replace = TRUE)

pvals1 <- numeric()
pvals2 <- numeric()
pvals3 <- numeric()
R2_means <- numeric()
R2_sds <- numeric()
pred_perf_R2s <- numeric()
x <- matrix(rnorm(n_samples*250), nrow = n_samples)
for (i in 1:500) {
  y <- rnorm(n_samples)
  fit <- cv.glmnet(x, y, nfolds = n_k_folds, alpha = 0.5, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
  lambda_seq <- fit$lambda
  cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
  fits <- list()
  fits[[1]] <- fit
  cvms <- matrix(nrow = 100, ncol = n_times)
  cvms[1:length(fit$cvm),1] <- fit$cvm
  #for (i in 2:(n_times)) {
  #  fit <- cv.glmnet(x, y, lambda = lambda_seq, nfolds = n_k_folds, alpha = 0.5, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
  #  fits[[i]] <- fit
  #  cvms[1:length(fit$cvm),i] <- fit$cvm
  #}
  avg_cvm <- rowMeans(cvms)
  best_lam_ind <- which.min(avg_cvm)
  best_lambda <- lambda_seq[best_lam_ind]
  performance <- evaluate_performance(x, y, fits[[1]], best_lam_ind, cv_fold_ids, n_k_folds)
  #pvals1 <- c(pvals1, performance$p_val)
  pvals2 <- c(pvals2, performance$pred_perf_pval)
  R2_means <- c(R2_means, performance$R2_mean)
  R2_sds <- c(R2_sds, performance$R2_sd)
  pred_perf_R2s <- c(pred_perf_R2s, performance$pred_perf_rsq)
}

#hist(pvals1, breaks = 20, main = 'P-val for F-Test - CV Elastic Net, Null')
hist(pvals2, breaks = 20, main = 'Predicted Performance P-val - CV Elastic Net, Null')


# New way ----
generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

library(doParallel)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores, type="PSOCK")
registerDoParallel(cl)
n_samples <- 150
n_train_test_folds <- 5
n_cv_folds <- 10

simulation_null <- function(n_samples, n_train_test_folds, n_cv_folds) {
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  x <- matrix(rnorm(n_samples*250), nrow = n_samples)
  y <- rnorm(n_samples)
  y_fold_pred <- rep(0, n_samples)
  R2_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    x_train <- x[which(train_test_fold_ids != test_fold), ]
    y_train <- y[which(train_test_fold_ids != test_fold)]
    x_test <- x[which(train_test_fold_ids == test_fold), ]
    y_test <- y[which(train_test_fold_ids == test_fold)]
    cv_fold_ids <- generate_fold_ids(length(y_train), n_cv_folds)
    fit <- cv.glmnet(x_train, y_train, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, parallel = TRUE)
    y_pred <- predict(fit, x_test, s = 'lambda.min')
    tss <- sum((y_test - mean(y_test))**2)
    rss <- sum((y_test - y_pred)**2)
    R2 <- 1 - rss/tss
    pval <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
    R2_folds[test_fold] <- R2
    pval_folds[test_fold] <- pval
    y_fold_pred[which(train_test_fold_ids == test_fold)] <- y_pred
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  correlation_test <- cor.test(y, y_fold_pred)
  R2_one_shot <- (correlation_test$estimate)**2
  pval_one_shot <- correlation_test$p.value
  c(R2_avg, R2_sd, pval_est, R2_one_shot, pval_one_shot)
}

simulation_signal <- function(n_samples, n_train_test_folds, n_cv_folds) {
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  x <- matrix(rnorm(n_samples*250), nrow = n_samples)
  signal <- FALSE
  if (runif(1) < 0.5) {
    y <- sqrt(0.2)*x[,1] + sqrt(0.3)*x[,2] + rnorm(n_samples, sd = sqrt(0.5))
    signal <- TRUE
  } else {
    y <- rnorm(n_samples)
  }
  y_fold_pred <- rep(0, n_samples)
  R2_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    x_train <- x[which(train_test_fold_ids != test_fold), ]
    y_train <- y[which(train_test_fold_ids != test_fold)]
    x_test <- x[which(train_test_fold_ids == test_fold), ]
    y_test <- y[which(train_test_fold_ids == test_fold)]
    cv_fold_ids <- generate_fold_ids(length(y_train), n_cv_folds)
    fit <- cv.glmnet(x_train, y_train, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, parallel = TRUE)
    y_pred <- predict(fit, x_test, s = 'lambda.min')
    tss <- sum((y_test - mean(y_test))**2)
    rss <- sum((y_test - y_pred)**2)
    R2 <- 1 - rss/tss
    pval <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
    R2_folds[test_fold] <- R2
    pval_folds[test_fold] <- pval
    y_fold_pred[which(train_test_fold_ids == test_fold)] <- y_pred
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  correlation_test <- cor.test(y, y_fold_pred)
  R2_one_shot <- (correlation_test$estimate)**2
  pval_one_shot <- correlation_test$p.value
  
  c(R2_avg, R2_sd, pval_est, R2_one_shot, pval_one_shot, signal)
}

simulation_signal2 <- function(n_samples, n_train_test_folds, n_cv_folds) {
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  x <- matrix(rnorm(n_samples*250), nrow = n_samples)
  signal <- FALSE
  if (runif(1) < .5) {
    y <- sqrt(0.2)*x[,1] + sqrt(0.3)*x[,2] + rnorm(n_samples, sd = sqrt(0.5))
    signal <- TRUE
  } else {
    y <- rnorm(n_samples)
  }
  y_fold_pred <- rep(0, n_samples)
  R2_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    x_train <- x[which(train_test_fold_ids != test_fold), ]
    y_train <- y[which(train_test_fold_ids != test_fold)]
    x_test <- x[which(train_test_fold_ids == test_fold), ]
    y_test <- y[which(train_test_fold_ids == test_fold)]
    cv_fold_ids <- generate_fold_ids(length(y_train), n_cv_folds)
    fit <- cv.glmnet(x_train, y_train, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, parallel = TRUE, keep = TRUE)
    if (test_fold == 1) {
      cross_val_summary <- summary(lm(y_train ~ fit$fit.preval[,which.min(fit$cvm)]*as.factor(cv_fold_ids)))
      cross_val_pval <- cross_val_summary$coef[2,4]
      cross_val_R2 <- cross_val_summary$r.squared
    }
    y_pred <- predict(fit, x_test, s = 'lambda.min')
    tss <- sum((y_test - mean(y_test))**2)
    rss <- sum((y_test - y_pred)**2)
    R2 <- 1 - rss/tss
    pval <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
    R2_folds[test_fold] <- R2
    pval_folds[test_fold] <- pval
    y_fold_pred[which(train_test_fold_ids == test_fold)] <- y_pred
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  
  one_shot_summary <- summary(lm(y ~ y_fold_pred*as.factor(train_test_fold_ids)))
  R2_one_shot <- one_shot_summary$r.squared
  pval_one_shot <- one_shot_summary$coef[2,4]
  c(R2_avg, R2_sd, pval_est, R2_one_shot, pval_one_shot, cross_val_R2, cross_val_pval, signal)
}


null_results <- matrix(ncol = 5, nrow = 1000)
for (i in 1:1000) {
  null_results[i,] <- simulation_null(n_samples, n_train_test_folds, n_cv_folds)
}

signal_results <- matrix(ncol = 6, nrow = 1000)
for (i in 1:1000) {
  signal_results[i,] <- simulation_signal(n_samples, n_train_test_folds, n_cv_folds)
}

cl <- makeCluster(n_cores, type="PSOCK")
registerDoParallel(n_cores-1)
signal_results2 <- matrix(ncol = 8, nrow = 1000)
for (i in 1:1000) {
  signal_results2[i,] <- simulation_signal2(n_samples, n_train_test_folds, n_cv_folds)
}


stopImplicitCluster()


