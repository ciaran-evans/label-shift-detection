library(mvtnorm)
library(dplyr)


# check this
args <- commandArgs(trailingOnly = TRUE)
ntrain_classifier = as.numeric(args[1])
batch_num = as.numeric(args[2])


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

qda_predict <- function(x, mu0_est, mu1_est, sig0_est, sig1_est, p1_est){
  return(p1_est*dmvnorm(x, mu1_est, sig1_est)/(p1_est*dmvnorm(x, mu1_est, sig1_est) +
                                                 (1 - p1_est)*dmvnorm(x, mu0_est, sig0_est)))
}


cusum_reduce <- function(w, l){
  return(max(0, w) + l)
}

# data_func generates data, 
# llr_func calculates log likelihood ratio
estimate_detection_time <- function(data_func, llr_func, nreps,
                                    cs_cutoff_list){
  stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
  for(rep in 1:nreps){
    obs_data <- data_func()
    llrs <- llr_func(obs_data)
    cs_stat <- Reduce(cusum_reduce, llrs, accumulate=T)
    for(cc in 1:length(cs_cutoff_list)){
      stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
    }
    
    print(rep)
  }
  mean_stop_time <- c()
  inf_stop_time <- c()
  for(cc in 1:length(cs_cutoff_list)){
    cc_stopping_times <- stopping_times[,cc]
    mean_stop_time[cc] <- mean(cc_stopping_times[cc_stopping_times < Inf])
    inf_stop_time[cc] <- sum(cc_stopping_times == Inf)
  }
  detection_output <- list(mean_stop_time = mean_stop_time,
                           inf_stop_time = inf_stop_time)
  return(detection_output)
}



d=150
pi_inf = 0.4
pi_0 = 0.7
mu0 = rep(0, d)
mu1 = rep(1, d)
sig0 = diag(d)
sig1 = matrix(0.1, d, d)
diag(sig1) = 2
nreps = 200
nobs_arl = 4000
nobs_dd = 800


if(ntrain_classifier == 800){
  lda_cutoff <- 2.9
  qda_cutoff <- 4.8
  seed_num <- batch_num
} else if(ntrain_classifier == 1000){
  lda_cutoff <- 3
  qda_cutoff <- 4
  seed_num <- 200 + batch_num 
} else if(ntrain_classifier == 1500){
  lda_cutoff <- 3
  qda_cutoff <- 3.3
  seed_num <- 400 + batch_num
} else if(ntrain_classifier == 2000){
  lda_cutoff <- 3
  qda_cutoff <- 3.2
  seed_num <- 600 + batch_num
} else if(ntrain_classifier == 3000){
  lda_cutoff <- 3
  qda_cutoff <- 3.2
  seed_num <- 800 + batch_num
} else if(ntrain_classifier == 5000){
  lda_cutoff <- 3
  qda_cutoff <- 3.2
  seed_num <- 1000 + batch_num
} else {
  lda_cutoff <- 2.8
  qda_cutoff <- 5
  seed_num <- 1200 + batch_num
}


lda_cutoff <- seq(lda_cutoff - 0.5, lda_cutoff + 0.5, length.out = 100)
qda_cutoff <- seq(qda_cutoff - 0.5, qda_cutoff + 0.5, length.out = 100)

set.seed(seed_num)

# Classifier training data
classifier_training_data <- generate_raw_data(d, ntrain_classifier, 
                                              pi_inf, mu0, mu1, sig0, sig1)
y_train_classifier <- classifier_training_data$ys
x_train_classifier <- classifier_training_data$xs


# LDA estimation

mu0_est <- colMeans(x_train_classifier[y_train_classifier == 0,])
mu1_est <- colMeans(x_train_classifier[y_train_classifier == 1,])

sig_est <- ((sum(y_train_classifier == 1) - 1)*cov(x_train_classifier[y_train_classifier == 1,]) + 
              (sum(y_train_classifier == 0) - 1)*cov(x_train_classifier[y_train_classifier == 0,]))/(ntrain_classifier - 1)

# QDA estimation
# things are the same except we estimate separate covariance matrices
sig0_est <- cov(x_train_classifier[y_train_classifier == 0,])
sig1_est <- cov(x_train_classifier[y_train_classifier == 1,])



# lda

llr_func <- function(x){
  lr <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf)) * lda_predict(x, mu0_est, mu1_est, sig_est, pi_inf) + (1 - pi_0)/(1 - pi_inf)
  return(log(lr))
}

# arl
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_inf, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, lda_cutoff)
lda_arl_means <- results$mean_stop_time
lda_arl_inf <- results$inf_stop_time

# detection delay
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_0, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, lda_cutoff)
lda_dd_means <- results$mean_stop_time
lda_dd_inf <- results$inf_stop_time




# qda

llr_func <- function(x){
  lr <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf)) * qda_predict(x, mu0_est, mu1_est, sig0_est, sig1_est, pi_inf) + (1 - pi_0)/(1 - pi_inf)
  return(log(lr))
}

# arl
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_inf, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, qda_cutoff)
qda_arl_means <- results$mean_stop_time
qda_arl_inf <- results$inf_stop_time

# detection delay
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_0, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, qda_cutoff)
qda_dd_means <- results$mean_stop_time
qda_dd_inf <- results$inf_stop_time

output <- data.frame(arl = c(lda_arl_means, qda_arl_means),
                     dd = c(lda_dd_means, qda_dd_means),
                     arl_inf = c(lda_arl_inf, qda_arl_inf),
                     dd_inf = c(lda_dd_inf, qda_dd_inf),
                     type = rep(c("KDE LDA", "KDE QDA"), each=100),
                     batch_num = batch_num,
                     ntrain_classifier = ntrain_classifier,
                     thresh = c(lda_cutoff, qda_cutoff))

write.table(output, file=paste("../simulation_output/lda_qda_comparison_classifier_", 
                               ntrain_classifier, "_",
                               batch_num,
                               ".txt", sep=""),
            sep=" ", row.names = F, col.names = F)

