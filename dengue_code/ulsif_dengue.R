library(readr)
library(dplyr)
library(tidyr)
library(mgcv)
library(matrixStats)
library(gridExtra)
library(densratio)

dengue_original <- read_csv("../dengue_original.csv")

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

set.seed(3)
train_inds <- sample(1:nrow(dengue), 1000, replace = F)
test_inds <- setdiff(1:nrow(dengue), train_inds)
dengue_train <- dengue[train_inds,]
dengue_test <- dengue[test_inds,]


train_data <- dengue_train %>%
  select(Age, WBC)

nsim = 50
cs_cutoff_list <- seq(0.005, 2, 0.005)
mean_arls <- matrix(nrow=nsim, ncol=length(cs_cutoff_list))
mean_dds <- matrix(nrow=nsim, ncol=length(cs_cutoff_list))
for(sim in 1:nsim){
  sim_test_data <- train_data[dengue_train$Dengue == 1,] %>%
    sample_n(680, replace=T) %>%
    rbind(train_data[dengue_train$Dengue == 0,] %>%
            sample_n(320, replace=T)) 
  
  test_data <- dengue_test %>%
    select(Age, WBC)
  
  ratio_est <- densratio(test_data, train_data, "uLSIF",
                         alpha = 0)
  
  test_lrs <- ratio_est$compute_density_ratio(test_data)
  
  nreps <- 500
  arl_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
  dd_stopping_times <- matrix(nrow = nreps, ncol = length(cs_cutoff_list))
  nobs_arl = 20000
  nobs_dd = 20000
  
  class_0_preds <- test_lrs[dengue_test$Dengue == 0]
  class_1_preds <- test_lrs[dengue_test$Dengue == 1]
  
  pi_inf = mean(dengue_train$Dengue)
  pi_0 = 0.68
  
  cusum_reduce <- function(w, l){
    return(max(0, w) + l)
  }
  
  for(rep in 1:nreps){
    sample_ys <- rbinom(nobs_arl, 1, pi_inf)
    sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
      sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
    lrs <- sample_obs
    cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
    for(cc in 1:length(cs_cutoff_list)){
      arl_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
    }
    
    sample_ys <- rbinom(nobs_dd, 1, pi_0)
    sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
      sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
    lrs <- sample_obs
    cs_stat <- Reduce(cusum_reduce, log(lrs), accumulate=T)
    for(cc in 1:length(cs_cutoff_list)){
      dd_stopping_times[rep, cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
    }
    
    print(rep)
  }
  
  for(cc in 1:length(cs_cutoff_list)){
    cc_stopping_times <- arl_stopping_times[,cc]
    mean_arls[sim, cc] <- mean(cc_stopping_times[cc_stopping_times < Inf],
                          na.rm=T)
  }
  
  for(cc in 1:length(cs_cutoff_list)){
    cc_stopping_times <- dd_stopping_times[,cc]
    mean_dds[sim, cc] <- mean(cc_stopping_times[cc_stopping_times < Inf],
                         na.rm=T)
  }
}


write.table(mean_arls, 
            file="../dengue_output/ulsif_dengue_arls.txt",
            sep = " ", row.names = F, col.names = F)

write.table(mean_dds, 
            file="../dengue_output/ulsif_dengue_dds.txt",
            sep = " ", row.names = F, col.names = F)
