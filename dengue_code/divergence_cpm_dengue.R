library(readr)
library(dplyr)
library(tidyr)
library(mgcv)
library(matrixStats)
library(gridExtra)
library(cpm)

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
  select(Vomiting, Skin, BMI, Age, Temperature,
         WBC, HCT, PLT) %>%
  mutate(Vomiting = ifelse(Vomiting == "yes", 1, 0),
         Skin = ifelse(Skin == "yes", 1, 0))

test_data <- dengue_test %>%
  select(Vomiting, Skin, BMI, Age, Temperature,
         WBC, HCT, PLT) %>%
  mutate(Vomiting = ifelse(Vomiting == "yes", 1, 0),
         Skin = ifelse(Skin == "yes", 1, 0))

train_xbar <- colMeans(train_data)

calc_similarity <- function(train_xbar, new_obs){
  return(sum(train_xbar * new_obs)/(sqrt(sum(train_xbar^2)) + sqrt(sum(new_obs^2))))
}

apply_func <- function(r){
  calc_similarity(train_xbar, r)
}

test_dists <- apply(test_data, 1,
                    apply_func)

class_0_preds <- c(test_dists)[dengue_test$Dengue == 0]
class_1_preds <- c(test_dists)[dengue_test$Dengue == 1]

nreps <- 10000
nobs_dd = 1000
nobs_train = 1000

pi_inf = mean(dengue_train$Dengue)
pi_0 = 0.68

cpm_means <- c()
cpm_sds <- c()
for(arl in c(370, 500, 700, 1000, 2000)){
  detection_times <- c()
  for(rep in 1:nreps){
    sample_ys <- c(rbinom(nobs_train, 1, pi_inf), rbinom(nobs_dd, 1, pi_0))
    sample_obs <- sample(class_1_preds, nobs_train + 
                           nobs_dd, replace=T)*sample_ys +
      sample(class_0_preds, nobs_train + 
               nobs_dd, replace=T)*(1 - sample_ys)
    
    cpm_stop <- detectChangePoint(sample_obs,
                                  "Cramer-von-Mises",
                                  arl, nobs_train)
    detection_times[rep] <- cpm_stop$detectionTime - nobs_train
    print(paste(arl, rep))
  }
  
  cpm_means <- c(cpm_means, mean(detection_times[detection_times >= 0]))
  cpm_sds <- c(cpm_sds, 
               sd(detection_times[detection_times >= 0])/sqrt(nreps))
}

write.table(cbind(cpm_means, cpm_sds), 
            file="../dengue_output/divergence_cpm_dengue_results.txt",
            sep = " ", row.names = F, col.names = F)
