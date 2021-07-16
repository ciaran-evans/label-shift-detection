library(mvtnorm)
library(dplyr)
library(densratio)


# check this
args <- commandArgs(trailingOnly = TRUE)
ntrain_classifier = as.numeric(args[1])
sigma_ind = as.numeric(args[2])
batch_num = as.numeric(args[3])


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

true_predict <- function(x, mu0, mu1, sig0, sig1, p1){
  return(p1*dmvnorm(x, mu1, sig1)/(p1*dmvnorm(x, mu1, sig1) +
                                                 (1 - p1)*dmvnorm(x, mu0, sig0)))
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
  }
  mean_stop_time <- c()
  inf_stop_time <- c()
  sd_stop_time <- c()
  for(cc in 1:length(cs_cutoff_list)){
    cc_stopping_times <- stopping_times[,cc]
    mean_stop_time[cc] <- mean(cc_stopping_times[cc_stopping_times < Inf])
    inf_stop_time[cc] <- sum(cc_stopping_times == Inf)
    sd_stop_time[cc] <- sd(cc_stopping_times[cc_stopping_times < Inf])
  }
  detection_output <- list(mean_stop_time = mean_stop_time,
                           inf_stop_time = inf_stop_time,
                           sd_stop_time = sd_stop_time)
  return(detection_output)
}



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

nreps = 1000
nobs_arl = 10000
nobs_dd = 2000

lda_cutoff <- seq(1.5, 5, length.out = 100)
ulsif_cutoff <- seq(1.5, 5, length.out = 100)
set.seed(batch_num)

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

# lda lr func
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
lda_arl_sd <- results$sd_stop_time

# detection delay
data_func <- function(){
  return(generate_raw_data(d, nobs_dd, pi_0, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, lda_cutoff)
lda_dd_means <- results$mean_stop_time
lda_dd_inf <- results$inf_stop_time
lda_dd_sd <- results$sd_stop_time


## ulsif
post_train_sample <- x_train_classifier[y_train_classifier == 0,] %>%
  as.data.frame() %>%
  sample_n(ntrain_classifier*(1 - pi_0), replace=T) %>%
  rbind(x_train_classifier[y_train_classifier == 1,] %>%
          as.data.frame() %>%
          sample_n(ntrain_classifier*pi_0, replace=T)) %>%
  as.matrix()

# post_train_sample <- generate_raw_data(d, ntrain_classifier,
#                                        pi_0, mu0, mu1, sig0, sig1)$xs

ratio_est <- densratio(post_train_sample, x_train_classifier, 
                       alpha = 0)


# ulsif lr func
llr_func <- function(x){
  return(log(ratio_est$compute_density_ratio(x)))
}

# arl
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_inf, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, ulsif_cutoff)
ulsif_arl_means <- results$mean_stop_time
ulsif_arl_inf <- results$inf_stop_time
ulsif_arl_sd <- results$sd_stop_time

# dd
data_func <- function(){
  return(generate_raw_data(d, nobs_dd, pi_0, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, ulsif_cutoff)
ulsif_dd_means <- results$mean_stop_time
ulsif_dd_inf <- results$inf_stop_time
ulsif_dd_sd <- results$sd_stop_time


output <- data.frame(arl = c(lda_arl_means, ulsif_arl_means),
                     dd = c(lda_dd_means, ulsif_dd_means),
                     arl_inf = c(lda_arl_inf, ulsif_arl_inf),
                     dd_inf = c(lda_dd_inf, ulsif_dd_inf),
                     arl_sd = c(lda_arl_sd, ulsif_arl_sd),
                     dd_sd = c(lda_dd_sd, ulsif_dd_sd),
                     type = rep(c("LDA", "uLSIF"), each=100),
                     batch_num = batch_num,
                     ntrain_classifier = ntrain_classifier,
                     thresh = c(lda_cutoff, ulsif_cutoff))

write.table(output, file=paste("../simulation_output/lda_ulsif_comparison_", 
                               ntrain_classifier, "_",
                               sigma_ind, "_",
                               batch_num,
                               ".txt", sep=""),
            sep=" ", row.names = F, col.names = F)
