library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm)
library(gridExtra)


#############################################################
######### Figure 2 (Example 3)                      #########
######### (comparison between LDA and QDA in high   ######### 
#########  dimensions, as sample size changes)      #########
#############################################################


# function for generating data from a mixture of multivariate normals
generate_raw_data <- function(d, n_obs, p1, mu0, mu1, sig0, sig1){
  ys <- rbinom(n_obs, 1, p1)
  xs <- rmvnorm(n_obs, mean=mu1, sigma=sig1)*matrix(rep(ys, d), ncol=d) +
    rmvnorm(n_obs, mean=mu0, sigma=sig0)*matrix(rep(1-ys, d), ncol=d)
  new_data <- list(ys = ys, xs=xs)
  return(new_data)
}

# likelihood ratio estimate under label shift assumption, 
# using LDA classifier with estimated means and covariance matrix
# (see equations (5) and (12) in the paper)
lda_lr <- function(x, mu0_est, mu1_est, sig_est){
  return((pi_0*dmvnorm(x, mu1_est, sig_est) + 
            (1 - pi_0)*dmvnorm(x, mu0_est, sig_est))/(pi_inf*dmvnorm(x, mu1_est, sig_est) + 
                                                        (1 - pi_inf)*dmvnorm(x, mu0_est, sig_est)))
}

# likelihood ratio estimate under label shift assumption, 
# using QDA classifier with estimated means and covariance matrices
# (see equation (5) in the paper)
qda_lr <- function(x, mu0_est, mu1_est, sig0_est, sig1_est){
  return((pi_0*dmvnorm(x, mu1_est, sig1_est) + 
            (1 - pi_0)*dmvnorm(x, mu0_est, sig0_est))/(pi_inf*dmvnorm(x, mu1_est, sig1_est) + 
                                                         (1 - pi_inf)*dmvnorm(x, mu0_est, sig0_est)))
}

d=150 # dimension of the data in Example 3
pi_inf = 0.4 # pre-change label distribution ( P(Y = 1) )
pi_0 = 0.7 # post-change label distribution ( P(Y = 1) )
mu0 = rep(0, d) # mean for class Y = 0
mu1 = rep(1, d) # mean for class Y = 1
sig0 = diag(d) # covariance matrix for class Y = 0
sig1 = matrix(0.1, d, d) # covariance matrix for class Y = 1
diag(sig1) = 2

# L1 distances for the LDA and QDA likelihood ratio estimates
lda_l1_dist <- c()
qda_l1_dist <- c()

set.seed(47) # set a seed for reproducibility

# We consider training sample sizes from 800 to 5000
# For small samples, LDA should do better than QDA because of
# the bias--variance trade-off. For larger samples, QDA should do
# better because the QDA assumptions are satisfied
for(ntrain_classifier in c(600, 800, 1000, 1500, 2000, 3000, 5000)){
  
  # Remember that performance of the detection procedure depends on 
  # the likelihood ratio estimate, which depends on the training data
  # To compare average performance, we sample many training samples
  for(rep in 1:200){
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
    
    # Now generate a new sequence of data to evaluate performance
    eval_data <- generate_raw_data(d, 2000, 
                                   pi_inf, mu0, mu1, sig0, sig1)$xs
    
    # estimated likelihood ratios on new data, using LDA 
    lda_lrs <- lda_lr(eval_data, mu0_est, mu1_est, sig_est)
    
    # estimated likelihood ratios on new data, using QDA
    qda_lrs <- qda_lr(eval_data, mu0_est, mu1_est, 
                      sig0_est, sig1_est)
    
    # true likelihood ratios, using true means and covariances
    true_lrs <- qda_lr(eval_data, mu0, mu1, 
                       sig0, sig1)
    
    # compare estimated likelihood ratios to truth with 
    # expected l1 distance
    # Note that the l1 distance is with respect to the
    # distribution of X
    lda_l1_dist <- c(lda_l1_dist,
                     mean(abs(lda_lrs - true_lrs)))
    qda_l1_dist <- c(qda_l1_dist,
                     mean(abs(qda_lrs - true_lrs)))
    
    # print(paste(ntrain_classifier, rep))
  }
}


l1_dists <- data.frame(dist = c(lda_l1_dist, qda_l1_dist),
                       ntrain = rep(rep(c(600, 800, 1000, 1500, 2000, 3000, 5000),
                                        each = 200), 2),
                       type = rep(c("LDA", "QDA"), each=1400))

p1 <- l1_dists %>%
  group_by(ntrain, type) %>%
  summarize(mean_dist = mean(dist), 
            sd_dist = 2*sd(dist)/sqrt(200)) %>%
  ggplot(aes(x = ntrain, y = mean_dist, color = type)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_dist - sd_dist, 
                    ymax = mean_dist + sd_dist)) +
  theme_bw() +
  labs(x = "Size of Classifier Training Sample", 
       y = "Expected L1 Distance", 
       color = "")

# Now compare detection delay as a function of sample size
# These results were generated by 
# lda_qda_changepoint_classifier_comparison.R

lda_qda_results <- data.frame()
for(m in c(600, 800, 1000, 1500, 2000, 3000, 5000)){
  for(i in 1:200){
    temp <- read_delim(paste("simulation_output/lda_qda_comparison_classifier_",  
                             m, "_", i, ".txt", sep=""), 
                       " ", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)
    lda_qda_results <- lda_qda_results %>%
      rbind(temp)
    
    # print(paste(m, i))
  }
}

colnames(lda_qda_results) <- c("arl", "dd", "arl_inf", "dd_inf",
                               "type", "sim_num", "ntrain_classifier",
                               "thresh")

res <- lda_qda_results %>%
  group_by(type, ntrain_classifier, thresh) %>%
  summarize(arl_mean = mean(arl[arl < Inf], na.rm=T),
            dd_mean = mean(dd[dd < Inf], na.rm=T),
            arl_sd = sd(arl[arl < Inf], na.rm=T),
            dd_sd = sd(dd[dd < Inf], na.rm=T),
            count = sum(dd < Inf, na.rm=T)) %>%
  drop_na() %>%
  mutate(dist = abs(arl_mean - 180)) %>%
  group_by(type, ntrain_classifier) %>%
  filter(dist == min(dist)) %>%
  mutate(mean_sd = 2*dd_sd/sqrt(count),
         type = ifelse(type == "KDE LDA", "LDA", "QDA"))

p2 <- res %>%
  ggplot(aes(x = ntrain_classifier, y = dd_mean, color=type)) +
  geom_line() +
  geom_errorbar(aes(ymin = dd_mean - mean_sd, ymax = dd_mean + mean_sd)) +
  theme_bw() +
  labs(x = "Size of Classifier Training Sample",
       y = "Average Detection Delay",
       color = "")


pdf(file="lda_qda_classifier_comparison_plot.pdf", width = 12, height = 4)
grid.arrange(p1, p2, ncol=2)
dev.off()


#############################################################
######### Table 2 results                           #########
######### (results for simulations, scenario 1,     #########
#########  where label shift holds but sample size  #########
#########  and covariance matrix change)            #########
#############################################################

### Notation: 
###   sigma_ind = 1 --> covariance matrix = [[1, 0], [0, 1]]
###   sigma_ind = 2 --> covariance matrix = [[2, 0.1], [0.1, 2]]
###   sigma_ind = 3 --> covariance matrix = [[4, 0.5], [0.5, 4]]
###   ntrain = training sample size (either 200, 1000, or 5000)

### Results from classifier CUSUM (lda) and uLSIF CUSUM
### Note that each of these methods depends on a training sample, 
### which adds randomness to performance
### To assess performance, for each simulation setting we average
### the results across 50 different training samples
### These results were generated by classifier_vs_ulsif_comparison.R

# data frame to store results
lda_ulsif <- data.frame(arl = c(),
                        dd = c(),
                        arl_inf = c(),
                        dd_inf = c(),
                        arl_sd = c(),
                        dd_sd = c(),
                        type = c(),
                        batch_num = c(),
                        ntrain_classifier = c(),
                        thresh = c(),
                        sigma_ind = c())

# We consider each ntrain, sigma_ind pair
# For each ntrain, sigma_ind pair there are 50 files 
# (one for each randomly selected training sample)
for(ntrain in c(200, 1000, 5000)){
  for(sigma_ind in c(1,2,3)){
    for(i in 1:50){
      temp <- read_delim(paste("simulation_output/lda_ulsif_comparison_", 
                               ntrain, "_", sigma_ind, "_",
                               i, ".txt", sep=""), 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
      temp <- cbind(temp, rep(sigma_ind, nrow(temp)))
      lda_ulsif <- rbind(lda_ulsif, temp)
      
      #print(paste(ntrain, sigma_ind, i))
    }
  }
}

colnames(lda_ulsif) <- c("arl", "dd", "arl_inf", "dd_inf",
                         "arl_sd", "dd_sd", "type", "batch_num",
                         "ntrain_classifier", "thresh", "sigma_ind")

# For each training sample, we calculate the expected stopping time 
# at a sequence of detection thresholds
# The marginal expected stopping time is then the average across
# training samples
# To compare performance, we assess detection delay (dd) when the 
# expected time to false alarm (arl) is approximately 500
lda_ulsif <- lda_ulsif %>%
  group_by(sigma_ind, ntrain_classifier, type, thresh) %>%
  summarize(arl = mean(arl, na.rm=T),
            dd_mean = mean(dd, na.rm=T),
            dd_sd = sd(dd, na.rm=T)/sqrt(50)) %>%
  mutate(diff = abs(arl - 500)) %>%
  filter(diff == min(diff))

# classifier CUSUM results, Table 2
classifier_cusum_results <- lda_ulsif %>%
  ungroup() %>%
  filter(type == "LDA") %>%
  select(sigma_ind, ntrain_classifier, dd_mean, dd_sd)

# uLSIF CUSUM results, Table 2
ulsif_cusum_results <- lda_ulsif %>%
  ungroup() %>%
  filter(type == "uLSIF") %>%
  select(sigma_ind, ntrain_classifier, dd_mean, dd_sd)




### Results for CPM (classifier)
### This is performance of the CPM method for detecting
### label shift, using Cramer--von-Mises tests to test for
### a change in distribution in the sequence of classifier scores
### These results are generated by the file lda_cpm_detection.R
### (In lda_cpm_detection.R, the detection threshold is selected 
###  so that ARL is approximately 500)

### config_num:
###   1 -- 3: ntrain = 200, sigma_ind = 1-3
###   4 -- 6: ntrain = 1000, sigma_ind = 1-3
###   7 -- 9: ntrain = 5000, sigma_ind = 1-3

cpm_results <- data.frame(cpm_delay = c(),
                         cpm_detect = c(),
                         config_num = c())
for(config_num in 1:9){
  temp <- read_delim(paste("simulation_output/cpm_changepoint_results_config_",
                           config_num, ".txt", sep=""),
                     " ", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
  cpm_results <- cpm_results %>%
    rbind(temp)
}
colnames(cpm_results) <- c("cpm_delay", "cpm_detect", "config_num")

# Results for CPM (classifier), Table 2
cpm_results %>%
  group_by(config_num) %>%
  summarize(dd = mean(cpm_delay[cpm_detect]),
            sd = sd(cpm_delay[cpm_detect])/sqrt(sum(cpm_detect)))



### Results for CPM (divergence)
### This is performance of the CPM method for detecting
### label shift, using Cramer--von-Mises tests to test for
### a change in distribution in the sequence of cosine divergences
### These results are generated by the file divergence_cpm_detection.R
### (In divergence_cpm_detection.R, the detection threshold is selected 
###  so that ARL is approximately 500)

### config_num:
###   1 -- 3: ntrain = 200, sigma_ind = 1-3
###   4 -- 6: ntrain = 1000, sigma_ind = 1-3
###   7 -- 9: ntrain = 5000, sigma_ind = 1-3

cpm_divergence_results <- data.frame(cpm_delay = c(),
                          cpm_detect = c(),
                          config_num = c())
for(config_num in 1:9){
  temp <- read_delim(paste("simulation_output/cpm_changepoint_divergence_results_config_",
                           config_num, ".txt", sep=""),
                     " ", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
  cpm_divergence_results <- cpm_divergence_results %>%
    rbind(temp)
}
colnames(cpm_divergence_results) <- c("cpm_delay", "cpm_detect", "config_num")

# Results for CPM (divergence), Table 2
cpm_divergence_results %>%
  group_by(config_num) %>%
  summarize(dd = mean(cpm_delay[cpm_detect]),
            sd = sd(cpm_delay[cpm_detect])/sqrt(sum(cpm_detect)))


### Results for kNN changepoint detection 
### For the kNN procedure, a window of size 200 is used, so
### we only consider a training sample of size 200.
### Thus, there are only three different simulation configurations
### for the kNN results:
###   config_1 --> sigma_ind = 1
###   config_2 --> sigma_ind = 2
###   config_3 --> sigma_ind = 3
### These results are generated by knn_changepoint_config*.R
### (In knn_changepoint_config*.R, the detection threshold is 
###  set so that ARL is approximately 500)

knn_results <- read_delim("simulation_output/knn_changepoint_results_config_1.txt", 
                          " ", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE) %>%
  rbind(read_delim("simulation_output/knn_changepoint_results_config_2.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_changepoint_results_config_3.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE))

colnames(knn_results) <- c("dd", "config_num")

# kNN results for Table 2
# (Note that these are lower bounds, as discussed in the caption 
#  to Table 2)
knn_results %>%
  mutate(dd = ifelse(dd < Inf, dd, 300)) %>%
  group_by(config_num) %>%
  summarize(dd_mean = mean(dd),
            dd_sd = sd(dd)/sqrt(1000))



### Results for Optimal CUSUM
### Recall that the optimal CUSUM procedure uses the 
### true pre- and post-change distributions for the
### detection statistic
### The code to calculate the expected stopping times is 
### provided here

# function to generate data from a mixture of multivariate gaussians
generate_raw_data <- function(d, n_obs, p1, mu0, mu1, sig0, sig1){
  ys <- rbinom(n_obs, 1, p1)
  xs <- rmvnorm(n_obs, mean=mu1, sigma=sig1)*matrix(rep(ys, d), ncol=d) +
    rmvnorm(n_obs, mean=mu0, sigma=sig0)*matrix(rep(1-ys, d), ncol=d)
  new_data <- list(ys = ys, xs=xs)
  return(new_data)
}

# true P(Y = 1 | X)
true_predict <- function(x, mu0, mu1, sig0, sig1, p1){
  return(p1*dmvnorm(x, mu1, sig1)/(p1*dmvnorm(x, mu1, sig1) +
                                     (1 - p1)*dmvnorm(x, mu0, sig0)))
}

# function for calculating (log) cusum statistic from sequence of
# likelihood ratios
cusum_reduce <- function(w, l){
  return(max(0, w) + l)
}

# function for estimating the average stopping time, 
# given a data-generating function, a function to compute the
# log likelihood ratio, a number of repetitions, and
# a sequence of detection thresholds
# For each detection threshold, results are averaged across
# the repetitions of the simulation
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
    #print(rep)
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

d=2 # dimensionality of data
pi_inf = 0.4 # pre-change P(Y = 1)
pi_0 = 0.7 # post-change P(Y = 1)
mu0 = rep(0, d) # mean for class Y = 0
mu1 = rep(1.5, d) # mean for class Y = 1
sig0 = diag(d) # covariance for class Y = 0
# sig1 = covariance for class Y = 1
# Note that there is no reliance on training samples, 
# so we only consider the different covariance matrices
sigma_ind = 1 # repeat for 2, 3
if(sigma_ind == 1){
  sig1 = diag(d)
} else if(sigma_ind == 2){
  sig1 = matrix(0.1, d, d)
  diag(sig1) = 2
} else {
  sig1 = matrix(0.5, d, d)
  diag(sig1) = 4
}

# number of repetitions for approximating average stopping time
nreps = 5000 

# length of sequences for calculating ARL and detection delay
# (more computationally efficient to specify a large fixed length 
#  than to use a while loop)
nobs_arl = 10000
nobs_dd = 2000

# detection thresholds to consider
cs_cutoff_list <- seq(3, 4, 0.01)

# function for calculating true log likelihood ratio of a new observation
llr_func <- function(x){
  lr <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf)) * true_predict(x, mu0, mu1, sig0, sig1, pi_inf) + (1 - pi_0)/(1 - pi_inf)
  return(log(lr))
}

set.seed(47) # set seed for reproducibility
# arl
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_inf, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, cs_cutoff_list)

arl_means <- results$mean_stop_time
arl_inf <- results$inf_stop_time
arl_sd <- results$sd_stop_time

# dd
data_func <- function(){
  return(generate_raw_data(d, nobs_dd, pi_0, mu0, mu1, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, cs_cutoff_list)

dd_means <- results$mean_stop_time
dd_inf <- results$inf_stop_time
dd_sd <- results$sd_stop_time

# get detection delay when ARL is approximately 500
dd_means[which.min(abs(arl_means - 500))]




#############################################################
######### Table 3 results                           #########
######### (results for simulations, scenario 2,     #########
#########  where sample size is fixed, but the      #########
#########  covariance matrix changes and the        #########
#########  label shift assumption is violated)      #########
#############################################################

### Notation: 
###   sigma_ind = 1 --> covariance matrix = [[1, 0], [0, 1]]
###   sigma_ind = 2 --> covariance matrix = [[2, 0.1], [0.1, 2]]
###   sigma_ind = 3 --> covariance matrix = [[4, 0.5], [0.5, 4]]
###   mu_ind = 1 --> mu00 = [0.5, 0.5], mu01 = [1, 1]
###   mu_ind = 2 --> mu00 = [0.75, 0.75], mu01 = [0.75, 0.75]
###   mu_ind = 3 --> mu00 = [1, 1], mu01 = [0.5, 0.5]

### Results from classifier CUSUM (lda) and uLSIF CUSUM
### Note that each of these methods depends on a training sample, 
### which adds randomness to performance
### To assess performance, for each simulation setting we average
### the results across 50 different training samples
### These results were generated by classifier_vs_ulsif_break_label_shift.R

# data frame to store results
lda_ulsif_break <- data.frame(dd = c(),
                              dd_inf = c(),
                              dd_sd = c(),
                              type = c(),
                              batch_num = c(),
                              mu_ind = c(),
                              sigma_ind = c(),
                              thresh = c())

# We consider each mu_ind, sigma_ind pair
# For each mu_ind, sigma_ind pair there are 50 files 
# (one for each randomly selected training sample)
for(mu in c(1, 2, 3)){
  for(sigma_ind in c(1,2,3)){
    for(i in 1:50){
      temp <- read_delim(paste("simulation_output/lda_ulsif_ls_break_", 
                               mu, "_", sigma_ind, "_",
                               i, ".txt", sep=""), 
                         " ", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
      lda_ulsif_break <- rbind(lda_ulsif_break, temp)
      
      #print(paste(mu, sigma_ind, i))
    }
  }
}

colnames(lda_ulsif_break) <- c("dd", "dd_inf", "dd_sd", "type", "batch_num",
                               "mu_ind", "sigma_ind", "thresh")


# For each training sample, we calculate the expected stopping time 
# at a sequence of detection thresholds
# The marginal expected stopping time is then the average across
# training samples
# To compare performance, we assess detection delay (dd) when the 
# expected time to false alarm (arl) is approximately 500;
# we get the ARL from the lda_ulsif data we imported earlier

lda_ulsif_break <- lda_ulsif_break %>%
  group_by(sigma_ind, mu_ind, type, thresh) %>%
  summarize(dd_mean = mean(dd, na.rm=T),
            dd_sd = sd(dd, na.rm=T)/sqrt(50))


lda_ulsif_break <- lda_ulsif_break %>%
  inner_join(lda_ulsif %>%
               filter(ntrain_classifier == 1000) %>%
               select(sigma_ind, type, thresh, arl), 
             by=c("sigma_ind", "type", "thresh"))

# classifier CUSUM results, Table 3
lda_ulsif_break %>%
  ungroup() %>%
  filter(type == "LDA") %>%
  select(sigma_ind, mu_ind, dd_mean, dd_sd)

# uLSIF CUSUM results, Table 3
lda_ulsif_break %>%
  ungroup() %>%
  filter(type == "uLSIF") %>%
  select(sigma_ind, mu_ind, dd_mean, dd_sd)



### Results for CPM (classifier)
### This is performance of the CPM method for detecting
### label shift, using Cramer--von-Mises tests to test for
### a change in distribution in the sequence of classifier scores
### These results are generated by the file 
###  lda_cpm_detection_break_label_shift.R
### (In lda_cpm_detection_break_label_shift.R, the detection threshold 
### is selected so that ARL is approximately 500)

### config_num:
###   1 -- 3: mu_ind = 1, sigma_ind = 1-3
###   4 -- 6: mu_ind = 2, sigma_ind = 1-3
###   7 -- 9: mu_ind = 3, sigma_ind = 1-3

cpm_results <- data.frame(cpm_delay = c(),
                          cpm_detect = c(),
                          config_num = c())
for(config_num in 1:9){
  temp <- read_delim(paste("simulation_output/cpm_changepoint_lda_break_ls_config_",
                           config_num, ".txt", sep=""),
                     " ", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
  cpm_results <- cpm_results %>%
    rbind(temp)
}
colnames(cpm_results) <- c("cpm_delay", "cpm_detect", "config_num")

# Results for CPM (classifier), Table 3
cpm_results %>%
  group_by(config_num) %>%
  summarize(dd = mean(cpm_delay[cpm_detect]),
            sd = sd(cpm_delay[cpm_detect])/sqrt(sum(cpm_detect)))


### Results for CPM (divergence)
### This is performance of the CPM method for detecting
### label shift, using Cramer--von-Mises tests to test for
### a change in distribution in the sequence of cosine divergences
### These results are generated by the file 
###  divergence_cpm_detection_break_label_shift.R
### (In divergence_cpm_detection_break_label_shift.R, the detection 
###  threshold is selected so that ARL is approximately 500)

### config_num:
###   1 -- 3: mu_ind = 1, sigma_ind = 1-3
###   4 -- 6: mu_ind = 2, sigma_ind = 1-3
###   7 -- 9: mu_ind = 3, sigma_ind = 1-3
cpm_divergence_results <- data.frame(cpm_delay = c(),
                                     cpm_detect = c(),
                                     config_num = c())
for(config_num in 1:9){
  temp <- read_delim(paste("simulation_output/cpm_changepoint_divergence_break_ls_config_",
                           config_num, ".txt", sep=""),
                     " ", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
  cpm_divergence_results <- cpm_divergence_results %>%
    rbind(temp)
}
colnames(cpm_divergence_results) <- c("cpm_delay", "cpm_detect", "config_num")

# Results for CPM (divergence), Table 3
cpm_divergence_results %>%
  group_by(config_num) %>%
  summarize(dd = mean(cpm_delay[cpm_detect]),
            sd = sd(cpm_delay[cpm_detect])/sqrt(sum(cpm_detect)))




### Results for kNN changepoint detection 
### For the kNN procedure, a window of size 200 is used, so
### we only consider a training sample of size 200
### config_num:
###   1 -- 3: mu_ind = 1, sigma_ind = 1-3
###   4 -- 6: mu_ind = 2, sigma_ind = 1-3
###   7 -- 9: mu_ind = 3, sigma_ind = 1-3
### These results are generated by knn_changepoint_break_ls_config*.R
### (In knn_changepoint_break_ls_config*.R, the detection threshold is 
###  set so that ARL is approximately 500)

knn_results <- read_delim("simulation_output/knn_results_break_ls_config_1.txt", 
                          " ", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_2.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_3.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_4.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_5.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_6.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_7.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_8.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE)) %>%
  rbind(read_delim("simulation_output/knn_results_break_ls_config_9.txt", 
                   " ", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE))

colnames(knn_results) <- c("dd", "config_num")

knn_results %>%
  mutate(dd = ifelse(dd < Inf, dd, 300)) %>%
  group_by(config_num) %>%
  summarize(dd_mean = mean(dd),
            dd_sd = sd(dd)/sqrt(1000))





### Results for Optimal CUSUM
### Recall that the optimal CUSUM procedure uses the 
### true pre- and post-change distributions for the
### detection statistic
### The code to calculate the expected stopping times is 
### provided here

# function to generate data from a mixture of multivariate gaussians
generate_raw_data <- function(d, n_obs, p1, mu0, mu1, sig0, sig1){
  ys <- rbinom(n_obs, 1, p1)
  xs <- rmvnorm(n_obs, mean=mu1, sigma=sig1)*matrix(rep(ys, d), ncol=d) +
    rmvnorm(n_obs, mean=mu0, sigma=sig0)*matrix(rep(1-ys, d), ncol=d)
  new_data <- list(ys = ys, xs=xs)
  return(new_data)
}

# function for calculating (log) cusum statistic from sequence of
# likelihood ratios
cusum_reduce <- function(w, l){
  return(max(0, w) + l)
}

# function for estimating the average stopping time, 
# given a data-generating function, a function to compute the
# log likelihood ratio, a number of repetitions, and
# a sequence of detection thresholds
# For each detection threshold, results are averaged across
# the repetitions of the simulation
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
    
    #print(rep)
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

config_num <- 1 # repeat with 2 - 9
mu_ind <- rep(c(1, 2, 3), each=3)[config_num]
sigma_ind <- rep(c(1, 2, 3), 3)[config_num]

d=2  # dimensionality of data
pi_inf = 0.4 # pre-change P(Y = 1)
pi_0 = 0.7 # post-change P(Y = 1)
mu0_pre = rep(0, d) # pre-change mean of class Y = 0
mu1_pre = rep(1.5, d) # pre-change mean of class Y = 1
# mu0_post = post-change mean of class Y = 0
# mu1_post = post-change mean of class Y = 1
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

sig0 = diag(d) # covariance for class Y = 0
# sig1 = covariance for class Y = 1
if(sigma_ind == 1){
  sig1 = diag(d)
} else if(sigma_ind == 2){
  sig1 = matrix(0.1, d, d)
  diag(sig1) = 2
} else {
  sig1 = matrix(0.5, d, d)
  diag(sig1) = 4
}

# number of repetitions for approximating average stopping time
nreps = 5000

# length of sequences for calculating ARL and detection delay
# (more computationally efficient to specify a large fixed length 
#  than to use a while loop)
nobs_arl = 10000
nobs_dd = 2000

# detection thresholds to consider
cs_cutoff_list <- seq(2, 4.5, 0.05)

# true lr func
llr_func <- function(x){
  lr <- (pi_0*dmvnorm(x, mu1_post, sig1) + 
           (1-pi_0)*dmvnorm(x, mu0_post, sig0))/(
             pi_inf*dmvnorm(x, mu1_pre, sig1) +
               (1 - pi_inf)*dmvnorm(x, mu0_pre, sig0)
           )
  return(log(lr))
}

set.seed(47) # set seed for reproducibility
# arl
data_func <- function(){
  return(generate_raw_data(d, nobs_arl, pi_inf, mu0_pre, mu1_pre, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, cs_cutoff_list)

arl_means <- results$mean_stop_time
arl_inf <- results$inf_stop_time
arl_sd <- results$sd_stop_time

# dd
data_func <- function(){
  return(generate_raw_data(d, nobs_dd, pi_0, mu0_post, mu1_post, sig0, sig1)$xs)
}
results <- estimate_detection_time(data_func, llr_func, nreps, cs_cutoff_list)

dd_means <- results$mean_stop_time
dd_inf <- results$inf_stop_time
dd_sd <- results$sd_stop_time

# get detection delay when ARL is approximately 500
dd_means[which.min(abs(arl_means - 500))]

