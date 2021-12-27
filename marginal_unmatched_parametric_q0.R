library(hydroGOF)
library(tidyverse)
library(splitstackshape)
library(boot)

rm(list = ls())

PrW <- 0.5
AonY <- 0.5
superN <- 1000000
simN <- 100

WonA <- 2
WonY <- 2

n <- 5000

W <- rbinom(superN, 1, PrW) 
A <- rbinom(superN,1,plogis(-2.40 + WonA*W))
PrA <- mean(A)
Y <- rbinom(superN,1,plogis(-2.46 + AonY*A + WonY*W))
PrY <- mean(Y)
Y1 <- rbinom(superN,1,plogis(-2.46 + AonY + WonY*W))
Y0 <- rbinom(superN,1,plogis(-2.46 + WonY*W))
R1 <- mean(Y1)
R0 <- mean(Y0)
ATE <- R1 - R0

dat <- data.frame(W,A,Y)
id <- rownames(dat)
dat <- cbind(id = id, dat)

matched <- function(obsdat){
  
  obsdat <- dat[sample(x = superN, size = n, replace = F), ]
  
  ## Calculate the risk difference in the cohort
  
  Y0_cohort <-
    obsdat %>%
    filter(A == 0) 
  
  R0_cohort <- mean(pull(Y0_cohort,Y))
  
  Y1_cohort <-
    obsdat %>%
    filter(A == 1) 
  
  R1_cohort <- mean(pull(Y1_cohort,Y))
  
  RD_cohort <- R1_cohort - R0_cohort
  
  ## Select unmatched controls
  
  cases <- subset(obsdat, Y == 1, select = W:Y)
  controls <- subset(obsdat, Y == 0, select = W:Y)
  
  index <- sample(1:nrow(controls), nrow(cases),replace = FALSE)
  controls_selected <- controls[index, ]
  
  index_cases <- sample(1:nrow(cases), nrow(cases)*PrY,replace = FALSE)
  cases_adjusted <- cases[index_cases, ]
  
  index_controls <- sample(1:nrow(controls_selected), nrow(controls_selected)*(1 - PrY),replace = FALSE)
  controls_adjusted <- controls[index_controls, ]
  
  total <- bind_rows(cases_adjusted, controls_adjusted)
  
  Y0_step1 <-
    total %>%
    filter(A == 0) 
  
  R0_step1 <- mean(pull(Y0_step1,Y))
  
  Y1_step1 <-
    total %>%
    filter(A == 1) 
  
  R1_step1 <- mean(pull(Y1_step1,Y))
  
  RD_step1 <- R1_step1 - R0_step1
  
  ## G-computation
  
  # create a dataset with 3 copies of each subject
  # 1st copy: equal to original one
  total <-
    total %>%
    mutate(interv = -1) 
  
  # 2nd copy: treatment set to 0, outcome to missing
  interv0 <- total %>%
    mutate(interv = 0, A = 0, Y = NA)
  
  # 3rd copy: treatment set to 1, outcome to missing
  interv1 <- total %>%
    mutate(interv = 1, A = 1, Y = NA)
  
  onesample <- bind_rows(total, interv0, interv1) # combining datasets
  
  # linear model to estimate mean outcome conditional on treatment and confounders
  # parameters are estimated using original observations only (total_predicted)
  # parameter estimates are used to predict mean outcome for observations with 
  # exposure set to 0 (interv=0) and to 1 (interv=1)
  
  std <- glm(Y ~ A*W, data = onesample)
  summary(std)   
  onesample$predicted_Y <- predict(std, onesample)
  
  # estimate mean outcome in each of the groups interv=0, and interv=1
  # this mean outcome is a weighted average of the mean outcomes in each combination 
  # of values of treatment and confounders, that is, the standardized outcome
  R0_step2 <- mean(onesample[which(onesample$interv == 0),]$predicted_Y)
  R1_step2 <- mean(onesample[which(onesample$interv == 1),]$predicted_Y)
  
  RD_step2 <- R1_step2 - R0_step2
  
  diff_1 <- RD_cohort - RD_step1
  diff_2 <- ATE - RD_step2
  
  # function to calculate difference in means
  standardization <- function(data,indices) {
    # create a dataset with 3 copies of each subject
    d <- data[indices,] # 1st copy: equal to original one`
    d$interv <- -1
    d0 <- d # 2nd copy: treatment set to 0, outcome to missing
    d0$interv <- 0
    d0$A <- 0
    d0$Y <- NA
    d1 <- d # 3rd copy: treatment set to 1, outcome to missing
    d1$interv <- 1
    d1$A <- 1
    d1$Y <- NA
    d.onesample <- rbind(d, d0, d1) # combining datasets
    
    # linear model to estimate mean outcome conditional on treatment and confounders
    # parameters are estimated using original observations only (interv= -1)
    # parameter estimates are used to predict mean outcome for observations with set 
    # treatment (interv=0 and interv=1)
    fit <- glm(Y ~ A*W, 
               data = d.onesample)
    
    d.onesample$predicted_meanY <- predict(fit, d.onesample)
    
    # estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1
    return(c(mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
               mean(d.onesample$predicted_meanY[d.onesample$interv == 0])))
  }
  
  # bootstrap 
  results_bootstrap <- boot(data = total, statistic = standardization, R = 100)
  
  mean <- results_bootstrap$t0
  
  # generating confidence intervals through se
  se <- c(sd(results_bootstrap$t[,1]))
  ll_se <- mean - qnorm(0.975)*se
  ul_se <- mean + qnorm(0.975)*se
  
  # generating confidence intervals through quantile
  
  ll_q <- quantile(results_bootstrap$t[,1], probs = c(0.025)) 
  ul_q <- quantile(results_bootstrap$t[,1], probs = c(0.975))
  
  bootstrap <- data.frame(cbind(c("Treatment - No Treatment"), mean, se,
                                ll_se, ul_se, ll_q, ul_q)) 
  
  bootstrap <-
    bootstrap %>%
    filter(V1 == "Treatment - No Treatment") %>%
    mutate(ll_se = as.numeric(ll_se), ul_se = as.numeric(ul_se),
           ll_q = as.numeric(ll_q), ul_q = as.numeric(ul_q)) %>%
    mutate(coverage_se = if_else(ll_se <= ATE & ul_se >= ATE, 1, 0),
           coverage_q = if_else(ll_q <= ATE & ul_q >= ATE, 1, 0))
  
  results <-
    tibble(RD_cohort,RD_step1,RD_step2,diff_1,diff_2,bootstrap)
  
  return(results)
}

diff_list <- lapply(1:simN, function(x) matched(obsdat))
diff_results <- do.call(rbind, diff_list)

final_results <- tibble(
  ATE,
  RD_cohort = mean(pull(diff_results,RD_cohort)),
  RD_step1 = mean(pull(diff_results,RD_step1)),
  RD_step2 = mean(pull(diff_results,RD_step2)),
  diff_1 = mean(pull(diff_results,diff_1)),
  diff_2 = mean(pull(diff_results,diff_2)),
  ll_se = mean(pull(diff_results,ll_se)),
  ul_se = mean(pull(diff_results,ul_se)),
  coverage_se = mean(pull(diff_results,coverage_se)),
  ll_q = mean(pull(diff_results,ll_q)),
  ul_q = mean(pull(diff_results,ul_q)),
  coverage_q = mean(pull(diff_results,coverage_q))
)

