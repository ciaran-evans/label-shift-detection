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

calc_similarity <- function(train_xbar, new_obs){
  return(sum(train_xbar * new_obs)/(sqrt(sum(train_xbar^2)) + sqrt(sum(new_obs^2))))
}


for(config_num in 1:9){
  ntrain_classifier <- rep(c(200, 1000, 5000), each=3)[config_num]
  sigma_ind <- rep(c(1, 2, 3), 3)[config_num]
  
  d=2
  pi_inf = 0.4
  pi_0 = 0.7
  mu0 = rep(0, d)
  mu1 = rep(1.5, d)
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
                                                  pi_inf, mu0, mu1, sig0, sig1)
    x_train_classifier <- classifier_training_data$xs
    
    train_xbar <- colMeans(x_train_classifier)
    
    apply_func <- function(r){
      calc_similarity(train_xbar, r)
    }
    
    lda_preds_0 <- c(apply(x_train_classifier, 1,
                           apply_func),
                     apply(generate_raw_data(d, nobs_post_cpm,
                                             pi_0, mu0, mu1, sig0, sig1)$xs, 1,
                           apply_func))
    
    res <- detectChangePoint(lda_preds_0, "Cramer-von-Mises", 500, startup = ntrain_classifier)
    res$detectionTime - ntrain_classifier
    cpm_delay[i] <- res$detectionTime - ntrain_classifier
    cpm_detect[i] <- res$changeDetected
    print(paste(config_num, i))
  }
  
  output <- data.frame(cpm_delay = cpm_delay,
                       cpm_detect = cpm_detect,
                       config_num = config_num)
  
  write.table(output, file=paste("cpm_changepoint_divergence_results_config_", 
                                 config_num, ".txt", sep=""),
              sep=" ", row.names = F, col.names = F)
}

