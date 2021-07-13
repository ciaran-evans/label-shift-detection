library(gStream)
library(mvtnorm)

generate_raw_data <- function(d, n_obs, p1, mu0, mu1, sig0, sig1){
  ys <- rbinom(n_obs, 1, p1)
  xs <- rmvnorm(n_obs, mean=mu1, sigma=sig1)*matrix(rep(ys, d), ncol=d) +
    rmvnorm(n_obs, mean=mu0, sigma=sig0)*matrix(rep(1-ys, d), ncol=d)
  new_data <- list(ys = ys, xs=xs)
  return(new_data)
}

config_num <- 7
mu_ind <- rep(c(1, 2, 3), each=3)[config_num]
sigma_ind <- rep(c(1, 2, 3), 3)[config_num]

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

nobs_pre_nn = 200
nobs_post_nn = 300

nreps_nn <- 1000

detect_delays <- c()
for(rep in 1:nreps_nn){
  x_0 <- rbind(generate_raw_data(d, nobs_pre_nn, pi_inf, 
                                 mu0_pre, mu1_pre, sig0, sig1)$xs,
               generate_raw_data(d, nobs_post_nn, pi_0, 
                                 mu0_post, mu1_post, sig0, sig1)$xs)
  
  x0_dist <- as.matrix(dist(x_0))
  diag(x0_dist) = max(x0_dist)+100
  # the trick is to make the diagonal of the distance matrix really big
  r <- gstream(x0_dist, L = nobs_pre_nn, N0 = nobs_pre_nn, k = 5, statistics = 'w', asymp = F, skew.corr = T,
               ARL = 500)
  zstats <- r$scanZ$weighted
  thresh <- r$b$weighted
  detect_delays[rep] <- min(which(zstats > thresh))
  print(rep)
}

output <- data.frame(detect_delays = detect_delays,
                     config_num = config_num)

write.table(output, file=paste("knn_results_break_ls_config_", 
                               config_num, ".txt", sep=""),
            sep=" ", row.names = F, col.names = F)
