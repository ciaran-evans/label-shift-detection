library(mvtnorm)
library(cpm)


# functions we need

generate_raw_data <- function(d, n_obs, p1, mu0, mu1, sig0, sig1){
  ys <- rbinom(n_obs, 1, p1)
  xs <- rmvnorm(n_obs, mean=mu1, sigma=sig1)*matrix(rep(ys, d), ncol=d) +
    rmvnorm(n_obs, mean=mu0, sigma=sig0)*matrix(rep(1-ys, d), ncol=d)
  new_data <- list(ys = ys, xs=xs)
  return(new_data)
}

lda_predict <- function(x, mu0_est, mu1_est, sig_est, p1_est){
  return(p1_est*dmvnorm(x, mu1_est, sig_est)/(p1_est*dmvnorm(x, mu1_est, sig_est) +
                                                (1 - p1_est)*dmvnorm(x, mu0_est, sig_est)))
}


for(config_num in 1:9){
  mu_ind <- rep(c(1, 2, 3), each=3)[config_num]
  sigma_ind <- rep(c(1, 2, 3), 3)[config_num]
  
  ntrain_classifier = 1000
  
  d=2
  pi_inf = 0.4
  pi_0 = 0.7
  mu0_pre = rep(0, d)
  mu1_pre = rep(1.5, d)
  if(mu_ind == 1){
    mu0_post = rep(0.5, d)
    mu1_post = rep(1, d)
  } else if(mu_ind == 2){
    mu0_post = rep(0.75, d)
    mu1_post = rep(0.75, d)
  } else {
    mu0_post = rep(1, d)
    mu1_post = rep(0.5, d)
  }
  
  sig0 = diag(d)
  if(sigma_ind == 1){
    sig1 = diag(d)
  } else if(sigma_ind == 2){
    sig1 = matrix(0.1, d, d)
    diag(sig1) = 2
  } else {
    sig1 = matrix(0.5, d, d)
    diag(sig1) = 4
  }
  
  nobs_post_cpm = 1000
  nreps_cpm = 1000
  
  cpm_delay <- c()
  cpm_detect <- c()
  for(i in 1:nreps_cpm){
    classifier_training_data <- generate_raw_data(d, ntrain_classifier, 
                                                  pi_inf, mu0_pre, mu1_pre, sig0, sig1)
    y_train_classifier <- classifier_training_data$ys
    x_train_classifier <- classifier_training_data$xs
    
    # LDA estimation
    
    mu0_est <- colMeans(x_train_classifier[y_train_classifier == 0,])
    mu1_est <- colMeans(x_train_classifier[y_train_classifier == 1,])
    
    sig_est <- ((sum(y_train_classifier == 1) - 1)*cov(x_train_classifier[y_train_classifier == 1,]) + 
                  (sum(y_train_classifier == 0) - 1)*cov(x_train_classifier[y_train_classifier == 0,]))/(ntrain_classifier - 1)
    
    
    lda_preds_0 <- c(lda_predict(x_train_classifier, mu0_est, mu1_est,
                                 sig_est, pi_inf),
                     lda_predict(generate_raw_data(d, nobs_post_cpm,
                                                   pi_0, mu0_post, mu1_post, sig0, sig1)$xs,
                                 mu0_est, mu1_est,
                                 sig_est, pi_inf))
    
    res <- detectChangePoint(lda_preds_0, "Cramer-von-Mises", 500, startup = ntrain_classifier)
    res$detectionTime - ntrain_classifier
    cpm_delay[i] <- res$detectionTime - ntrain_classifier
    cpm_detect[i] <- res$changeDetected
    print(paste(config_num, i))
  }
  
  output <- data.frame(cpm_delay = cpm_delay,
                       cpm_detect = cpm_detect,
                       config_num = config_num)
  
  write.table(output, file=paste("cpm_changepoint_lda_break_ls_config_", 
                                 config_num, ".txt", sep=""),
              sep=" ", row.names = F, col.names = F)
}

