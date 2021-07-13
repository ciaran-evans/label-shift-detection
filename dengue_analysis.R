library(readr)
library(dplyr)
library(tidyr)
library(mgcv)
library(ROCR)
library(matrixStats)
library(gridExtra)
library(RColorBrewer)
library(scales)

### Import data
### This is the original data from Tuan et al (2015)
dengue_original <- read_csv("dengue_original.csv")

### Do some data cleaning
dengue <- dengue_original %>%
  select(SiteNo, Sex, Age, DayDisease, Vomiting, Abdo, Muco, Skin,
         Temp, BMI, Height, Weight, Flush, Rash, WBC, HCT, PLT, 
         NS1_TRIP, Lab_Confirmed_Dengue) %>%
  rename(SiteNumber = SiteNo, 
         DiseaseDay = DayDisease,
         Abdominal = Abdo,
         Mucosal = Muco,
         Temperature = Temp, 
         RapidTest = NS1_TRIP,
         Dengue = Lab_Confirmed_Dengue) %>%
  mutate(Vomiting = ifelse(Vomiting == 1, "yes", "no"),
         Abdominal = ifelse(Abdominal == 1, "yes", "no"),
         Mucosal = ifelse(Mucosal == 1, "yes", "no"),
         Skin = ifelse(Skin == 1, "yes", "no"),
         Flush = ifelse(Flush == 1, "yes", "no"),
         Rash = ifelse(Rash == 1, "yes", "no"),
         RapidTest = ifelse(RapidTest == 1, "positive", "negative"),
         Dengue = 2 - Dengue,
         Sex = ifelse(Sex == 1, "male", "female")) %>% 
  drop_na()

set.seed(3) # set a seed for reproducibility

### Divide data into training and test data
train_inds <- sample(1:nrow(dengue), 1000, replace = F)
test_inds <- setdiff(1:nrow(dengue), train_inds)
dengue_train <- dengue[train_inds,]
dengue_test <- dengue[test_inds,]

### Fit a logistic GAM with the training data
dengue_gam <- dengue_train %>%
  gam(Dengue ~ Vomiting + Skin + s(BMI) + s(Age) +
        s(Temperature) + s(WBC) + s(HCT) + s(PLT),
      data = ., family = binomial())

### Apply logistic GAM to test data to get
### predicted probabilities
predictions <- predict.gam(dengue_gam, newdata = dengue_test,
                           type = "response")

### ROC curve (part of Figure 3)
pred <- prediction(c(predictions), dengue_test$Dengue)
perf <- performance(pred,"tpr","fpr")

p1 <- data.frame(fpr = perf@x.values[[1]],
           tpr = perf@y.values[[1]]) %>%
  ggplot(aes(x = fpr, y = tpr)) +
  geom_line() +
  labs(x = "False positive rate", y = "True positive rate",
       title = "Dengue classifier ROC curve (AUC = 0.83)") +
  theme_bw()

### In Section 5, we compare binarized predictions with two 
### different cutoffs; 0.5, and the 'optimal' cutoff which maximizes 
### sensitivity + specificity

perf@alpha.values[[1]][which.max((1 - perf@x.values[[1]]) + perf@y.values[[1]])]

### so optimal cutoff = 0.334


### AUC of the ROC curve (0.83)
performance(pred, "auc")@y.values

### function for calculating (log) CUSUM detection statistic from
### sequence of (log) likelihood ratios
cusum_reduce <- function(w, l){
  return(max(0, w) + l)
}


#############################################################
######## Detection performance for an abrupt change #########
#############################################################

### First consider the detection procedure which uses the
### logistic GAM predicted probabilities

# We estimate operating characteristics by averaging stopping times 
# across many repetitions of the simulation. To look at the relationship
# between ARL and detection delay, we calculate stopping times for a range
# of detection thresholds
cs_cutoff_list <- seq(0.5, 6, 0.05)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_arl = 20000
nobs_dd = 1000

# new sequences of observations will be simulated by bootstrapping 
# the observed data
class_0_preds <- c(predictions)[dengue_test$Dengue == 0]
class_1_preds <- c(predictions)[dengue_test$Dengue == 1]

pi_inf = mean(dengue_train$Dengue) # pre-change P(Y = 1 | X)
pi_0 = 0.68 # post-change P(Y = 1 | X)

# repeat simulation many times to estimate operating characteristics
for(rep in 1:nreps){
  
  # simulate a new sequence of pre-change observations for 
  # time to false alarm
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
    sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  # simulate a new sequence of observations for detection delay
  sample_ys <- rbinom(nobs_dd, 1, pi_0)
  sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
    sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

# Results for classifier CUSUM (predicted probability)
prob_arls <- colMeans(arl_stopping_times, na.rm=T)
prob_dds <- colMeans(dd_stopping_times, na.rm=T)

#plot(prob_arls, prob_dds, type="l", ylim = c(0, 45))




### Next, repeat the process with the detection procedure that uses 
### binarized predictions instead of predicted probabilities 
### Here we use the optimal threshold (0.334) which maximizes 
### sensitivity + specificity

cs_cutoff_list <- seq(0.8, 8, 0.05)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_arl = 20000
nobs_dd = 1000

class_0_preds <- (c(predictions) > 0.334)[dengue_test$Dengue == 0]
class_1_preds <- (c(predictions) > 0.334)[dengue_test$Dengue == 1]

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68

for(rep in 1:nreps){
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
    sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  sample_ys <- rbinom(nobs_dd, 1, pi_0)
  sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
    sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

# Results for classifier CUSUM (binary, threshold = 0.334)
opt_thresh_arls <- colMeans(arl_stopping_times, na.rm=T)
opt_thresh_dds <- colMeans(dd_stopping_times, na.rm=T)

#points(opt_thresh_arls, opt_thresh_dds, type="l", col="red")




### Next, repeat the process with the detection procedure that uses 
### binarized predictions instead of predicted probabilities 
### Here we use a threshold of 0.5

cs_cutoff_list <- seq(0.5, 6, 0.05)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_arl = 20000
nobs_dd = 1000

class_0_preds <- (c(predictions) > 0.5)[dengue_test$Dengue == 0]
class_1_preds <- (c(predictions) > 0.5)[dengue_test$Dengue == 1]

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68

for(rep in 1:nreps){
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
    sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  sample_ys <- rbinom(nobs_dd, 1, pi_0)
  sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
    sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

# Results for classifier CUSUM (binary, threshold = 0.5)
other_thresh_arls <- colMeans(arl_stopping_times)
other_thresh_dds <- colMeans(dd_stopping_times)

points(other_thresh_arls, other_thresh_dds, type="l", col="blue")




### Now calculate operating characteristics for the optimal CUSUM 
### procedure. The optimal procedure uses the true dengue labels

cs_cutoff_list <- seq(0.5, 6, 0.05)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_arl = 20000
nobs_dd = 1000

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68

for(rep in 1:nreps){
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_ys +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  sample_ys <- rbinom(nobs_dd, 1, pi_0)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_ys +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

# Results for optimal CUSUM (true dengue status)
optimal_arls <- colMeans(arl_stopping_times,na.rm=T)
optimal_dds <- colMeans(dd_stopping_times,na.rm=T)

#points(optimal_arls, optimal_dds, type="l", col="green")




### We also compare the performance using the binary results from 
### the NS1 rapid antigen test. The rapid test does not give 
### perfect classification so should do worse than using true dengue 
### status, but should be better than using classifier predictions

cs_cutoff_list <- seq(0.5, 7, 0.05)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_arl = 20000
nobs_dd = 1000

rapid_test <- ifelse(dengue_test$RapidTest == "negative", 0, 1)
class_0_preds <- rapid_test[dengue_test$Dengue == 0]
class_1_preds <- rapid_test[dengue_test$Dengue == 1]

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68

for(rep in 1:nreps){
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
    sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  sample_ys <- rbinom(nobs_dd, 1, pi_0)
  sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
    sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
  lrs <- (pi_0/pi_inf - (1-pi_0)/(1 - pi_inf))*sample_obs +
    (1 - pi_0)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

# Results for rapid test CUSUM
rapid_arls <- colMeans(arl_stopping_times)
rapid_dds <- colMeans(dd_stopping_times)

#points(rapid_arls, rapid_dds, type="l", col="orange")




### We also compare performance to label shift detection with CPM
### We consider both the CPM statistic using cosine divergences, 
### and the CPM statistic using the classifier predicted probabilities

# these results generated by the file divergence_cpm_dengue.R
divergence_cpm_dengue_results <- read_delim("dengue_output/divergence_cpm_dengue_results.txt", 
                                            " ", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)

# these results generated by the file classifier_cpm_dengue.R
classifier_cpm_dengue_results <- read_delim("dengue_output/classifier_cpm_dengue_results.txt", 
                                            " ", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE)




### So far, we have assumed that pi0 (0.68) is known
### If pi0 is unknown, we can easily mix over a range of
### potential values for pi0
### These results were generated by the file 
### dengue_mixture_cusum_abrupt.R
### For computational efficiency, this is run in parallel
### across 500 iterations, with an output file for each iteration

mixture_arls <- matrix(nrow=500, ncol=171)
mixture_dds <- matrix(nrow=500, ncol=171)
for(i in 1:500){
  mixture_arls[i,] <- c(read_csv(paste("dengue_output/dengue_mixture_cusum_abrupt_arls_",
                                       i, ".txt", sep=""), 
           col_names = FALSE))$X1
  mixture_dds[i,] <- c(read_csv(paste("dengue_output/dengue_mixture_cusum_abrupt_dds_",
                                      i, ".txt", sep=""), 
                                col_names = FALSE))$X1
}

# Results for mixture CUSUM (predicted probability)
abrupt_mix_arls <- colMeans(mixture_arls)
abrupt_mix_dds <- colMeans(mixture_dds)




### Finally, we consider uLSIF, but performance is very poor
### uLSIF results were generated by ulsif_dengue.R

ulsif_dengue_arls <- read_delim("dengue_output/ulsif_dengue_arls.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE) %>%
  as.matrix()

ulsif_dengue_dds <- read_delim("dengue_output/ulsif_dengue_dds.txt", 
                                " ", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE) %>%
  as.matrix()

ulsif_arls <- colMeans(ulsif_dengue_arls, na.rm=T)
ulsif_dds <- colMeans(ulsif_dengue_dds)


### Plot results for abrupt change (part of Figure 3)

p2 <- data.frame(arl = c(prob_arls, opt_thresh_arls, 
                   other_thresh_arls, optimal_arls,
                   rapid_arls, abrupt_mix_arls,
                   c(370, 500, 700, 1000, 2000),
                   c(370, 500, 700, 1000, 2000)),
           dd = c(prob_dds, opt_thresh_dds, 
                  other_thresh_dds, optimal_dds, 
                  rapid_dds, abrupt_mix_dds,
                  divergence_cpm_dengue_results$X1,
                  classifier_cpm_dengue_results$X1),
           type = c(rep("Classifier CUSUM \n (predicted probability)", length(prob_arls)),
                    rep("Classifier CUSUM \n (binary, threshold = 0.334)",
                        length(opt_thresh_arls)),
                    rep("Classifier CUSUM \n (binary, threshold = 0.5)",
                        length(other_thresh_dds)),
                    rep("Optimal CUSUM \n (true dengue status)", length(optimal_arls)),
                    rep("Rapid test CUSUM", length(rapid_arls)),
                    rep("Mixture CUSUM \n (predicted probability)", length(abrupt_mix_arls)),
                    rep("CPM (divergence)", 5),
                    rep("CPM (classifier)", 5))) %>%
  mutate(type = factor(type,
                       levels = c("CPM (divergence)",
                                  "Classifier CUSUM \n (binary, threshold = 0.5)",
                                  "CPM (classifier)",
                                  "Classifier CUSUM \n (binary, threshold = 0.334)",
                                  "Classifier CUSUM \n (predicted probability)",
                                  "Mixture CUSUM \n (predicted probability)",
                                  "Rapid test CUSUM",
                                  "Optimal CUSUM \n (true dengue status)"))) %>%
  #filter(arl <= 150) %>%
  ggplot(aes(x = arl, y = dd, color=type)) +
  geom_smooth(se=F) +
  scale_x_continuous(limits=c(0, 2000)) +
  scale_color_manual(values = hue_pal()(8)) +
  labs(x = "Average time to false alarm", y = "Average detection delay",
       color = "Detection Procedure", 
       title = "Detection performance for an abrupt change") +
  theme_bw()






#############################################################
######## Detection performance for a gradual change #########
#############################################################


### First, we consider the mixture procedure, which mixes over 
### a range of possible values for pi0

# when no change occurs, the arls are the same as what we 
# already calculated above
gradual_mix_arls <- abrupt_mix_arls

# detection delays are calculated by the file dengue_mixture_cusum_gradual_v2.R
gradual_mix_dds <- read_delim("dengue_output/dengue_mixture_cusum_gradual_dds_v2.txt",
                              " ", escape_double = FALSE, 
                              col_names = FALSE)
gradual_mix_dds <- colMeans(gradual_mix_dds)



### In the gradual change case, we also consider the CPM procedure, 
### which performed well in the abrupt case, and is a natural choice 
### because it does not require knowledge of the post-change distribution

# these results are generated by classifier_cpm_dengue_gradual.R
gradual_change_cpm <- read_delim("dengue_output/classifier_cpm_dengue_gradual_results.txt", 
           " ", escape_double = FALSE, col_names = FALSE, 
           trim_ws = TRUE)





### Finally, we compare to the optimal CUSUM procedure for a gradual 
### change, which uses true dengue status and the true P(Y = 1 | X) 
### at each time point

cs_cutoff_list <- seq(0.5, 11, 0.01)
nreps <- 10000
arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
nobs_dd = 2000
nobs_arl = 30000

pi_inf = mean(dengue_train$Dengue) # pre-change P(Y = 1 | X)

# sequence of P(Y = 1 | X) after change occurs (they evolve over time)
prob_sequence <- 0.13 + 0.6*exp(100*(1:nobs_dd)/nobs_dd - 1)/(1 + exp(100*(1:nobs_dd)/nobs_dd - 1))

for(rep in 1:nreps){
  sample_ys <- rbinom(nobs_arl, 1, pi_inf)
  lrs <- (prob_sequence/pi_inf - (1-prob_sequence)/(1 - pi_inf))*sample_ys +
    (1 - prob_sequence)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  sample_ys <- rbinom(nobs_dd, 1, prob_sequence)
  lrs <- (prob_sequence/pi_inf - (1-prob_sequence)/(1 - pi_inf))*sample_ys +
    (1 - prob_sequence)/(1 - pi_inf)
  cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
  for(cc in 1:length(cs_cutoff_list)){
    dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
  }
  
  print(rep)
}

optimal_gradual_arls <- colMeans(arl_stopping_times, na.rm=T)
optimal_gradual_dds <- colMeans(dd_stopping_times, na.rm=T)


### Plot the results for a gradual change (part of Figure 3)

p3 <- data.frame(arl = c(gradual_mix_arls, optimal_gradual_arls,
                         c(370, 500, 700, 1000, 2000)),
           dd = c(gradual_mix_dds, optimal_gradual_dds,
                  gradual_change_cpm$X1),
           type = c(rep("Mixture CUSUM \n (predicted probability)", length(gradual_mix_arls)),
                    rep("Optimal CUSUM \n (true dengue status)", length(optimal_gradual_arls)),
                    rep("CPM (classifier)", 5))) %>%
  filter(arl <= 2000) %>%
  ggplot(aes(x = arl, y = dd, color = type)) +
  geom_line(lwd=1) +
  #scale_x_continuous(limits=c(0, 1500)) +
  labs(x = "Average time to false alarm", y = "Average detection delay",
       color = "Detection Procedure", 
       title = "Detection performance for a gradual change") +
  theme_bw() +
  scale_color_manual(values = hue_pal()(8)[c(3, 6, 8)])

layout_mat <- matrix(c(1, 1, 2, 2, 2,
                       1, 1, 2, 2, 2,
                       3, 3, 3, NA, NA,
                       3, 3, 3, NA, NA), ncol = 5, byrow = T)

pdf("dengue_detection_plots.pdf",
    width = 12, height=8)
grid.arrange(p1, p2, p3, layout_matrix = layout_mat)
dev.off()
