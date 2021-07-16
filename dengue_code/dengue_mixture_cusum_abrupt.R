library(readr)
library(dplyr)
library(tidyr)
library(mgcv)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
rep = as.numeric(args[1])

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

dengue_gam <- dengue_train %>%
  gam(Dengue ~ Vomiting + Skin + s(BMI) + s(Age) +
        s(Temperature) + s(WBC) + s(HCT) + s(PLT),
      data = ., family = binomial())

predictions <- predict.gam(dengue_gam, newdata = dengue_test,
                           type = "response")

window_size = 200
cs_cutoff_list <- exp(seq(0.5, 9, 0.05))
arl_stopping_times <- c()
dd_stopping_times <- c()
nobs_arl = 20000
nobs_dd = 1000


class_0_preds <- c(predictions)[dengue_test$Dengue == 0]
class_1_preds <- c(predictions)[dengue_test$Dengue == 1]

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68
num_pis <- 100

set.seed(rep)

sample_pis <- runif(num_pis, min=0.6, max = 0.8)
sample_ys <- rbinom(nobs_arl, 1, pi_inf)
sample_obs <- sample(class_1_preds, nobs_arl, replace=T)*sample_ys +
  sample(class_0_preds, nobs_arl, replace=T)*(1 - sample_ys)
lrs_mat <- matrix(nrow=num_pis, ncol=nobs_arl)
for(r in 1:num_pis){
  lrs_mat[r,] <- (sample_pis[r]/pi_inf - (1-sample_pis[r])/(1 - pi_inf))*sample_obs +
    (1 - sample_pis[r])/(1 - pi_inf)
}
cs_stat <- c()
for(t in 1:nobs_arl){
  window_vals <- c()
  for(k in max(t - window_size, 1):t){
    window_vals <- c(window_vals,
                     ifelse(k == t, mean(lrs_mat[,k]), 
                            mean(rowProds(lrs_mat[,k:t]))))
  }
  cs_stat[t] <- max(window_vals)
}
for(cc in 1:length(cs_cutoff_list)){
  arl_stopping_times[cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
}

sample_ys <- rbinom(nobs_dd, 1, pi_0)
sample_obs <- sample(class_1_preds, nobs_dd, replace=T)*sample_ys +
  sample(class_0_preds, nobs_dd, replace=T)*(1 - sample_ys)
lrs_mat <- matrix(nrow=num_pis, ncol=nobs_dd)
for(r in 1:num_pis){
  lrs_mat[r,] <- (sample_pis[r]/pi_inf - (1-sample_pis[r])/(1 - pi_inf))*sample_obs +
    (1 - sample_pis[r])/(1 - pi_inf)
}
cs_stat <- c()
for(t in 1:nobs_dd){
  window_vals <- c()
  for(k in max(t - window_size, 1):t){
    window_vals <- c(window_vals,
                     ifelse(k == t, mean(lrs_mat[,k]), 
                            mean(rowProds(lrs_mat[,k:t]))))
  }
  cs_stat[t] <- max(window_vals)
}
for(cc in 1:length(cs_cutoff_list)){
  dd_stopping_times[cc] <- min(which(cs_stat > cs_cutoff_list[cc]))
}

write.table(arl_stopping_times, file=
              paste("../dengue_output/dengue_mixture_cusum_abrupt_arls_",
                    rep, ".txt", sep=""),
      sep = " ", row.names = F, col.names = F)

write.table(dd_stopping_times, 
            paste("../dengue_output/dengue_mixture_cusum_abrupt_dds_",
                  rep, ".txt", sep=""),
            sep = " ", row.names = F, col.names = F)