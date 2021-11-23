library(hydroGOF)
library(tidyverse)
library(splitstackshape)
library(boot)

rm(list = ls())

PrW <- 0.5
AonY <- 0.5
superN <- 1000000
simN <- 10

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
  
  ## Select matched controls
  
  cases <- subset(obsdat, Y == 1)
  controls <- subset(obsdat, Y == 0)
  
  prop_W <- 
    controls %>%
    group_by(W) %>%
    summarise(prop = n()/nrow(controls))
  
  prop_W <- pull(prop_W, prop)
  prop_W0 <- prop_W[1]
  prop_W1 <- prop_W[2]
  
  cases_freq <- group_by(cases, W) %>% summarise(n = n_distinct(id))
  
  cs <- NA
  
  for (i in 1:nrow(cases_freq)) {
    pts <- controls[controls$W == cases_freq$W[i] ,]$id 
    recs <- sample(x = 1:length(pts), size = (cases_freq$n[i])) 
    cs <- c(cs, pts[recs]) 
  }
  
  cs <- cs[2:length(cs)] 
  cs <- as.data.frame(matrix(cs, ncol = 1))
  colnames(cs)[1] <- 'id'
  
  controls_selected <- left_join(x = cs,
                  y = controls, 
                  by = ('id' = 'id'))
  
  match <- rbind(cases,controls_selected)
  
  ## Create a pseudo cohort from the case-control sample
  
  con_W0_A <- 
    controls_selected %>%
    filter(W == 0)
  
  prop_W0_A <-
    con_W0_A %>%
    group_by(A) %>%
    summarise(prop = n()/nrow(con_W0_A))
  
  prop_W0_A <- pull(prop_W0_A, prop)
  prop_W0_A0 <- prop_W0_A[1]
  prop_W0_A1 <- prop_W0_A[2]
  
  con_W1_A <- 
    controls_selected %>%
    filter(W == 1)
  
  prop_W1_A <-
    con_W1_A %>%
    group_by(A) %>%
    summarise(prop = n()/nrow(con_W1_A))
  
  prop_W1_A <- pull(prop_W1_A, prop)
  prop_W1_A0 <- prop_W1_A[1]
  prop_W1_A1 <- prop_W1_A[2]
  
  type_1 <-
    tibble(W = 0, A = 0, Y = 0) %>%
    expandRows(type_1, count = prop_W0 * prop_W0_A0 * nrow(controls) , count.is.col = FALSE)
  
  type_2 <-
    tibble(W = 0, A = 1, Y = 0) %>%
    expandRows(type_2, count = prop_W0 * prop_W0_A1 * nrow(controls) , count.is.col = FALSE)
  
  type_3 <-
    tibble(W = 1, A = 0, Y = 0) %>%
    expandRows(type_3, count = prop_W1 * prop_W1_A0 * nrow(controls) , count.is.col = FALSE)
  
  type_4 <-
    tibble(W = 1, A = 1, Y = 0) %>%
    expandRows(type_4, count = prop_W1 * prop_W1_A1 * nrow(controls) , count.is.col = FALSE)
  
  controls_adjusted <-
    bind_rows(type_1, type_2, type_3, type_4)
  
  total <- bind_rows(cases, controls_adjusted)
  
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
    return(c(mean(d.onesample$predicted_meanY[d.onesample$interv == -1]),
             mean(d.onesample$predicted_meanY[d.onesample$interv == 0]),
             mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
             mean(d.onesample$predicted_meanY[d.onesample$interv == 1]) -
               mean(d.onesample$predicted_meanY[d.onesample$interv == 0])))
  }
  
  # bootstrap 
  results_bootstrap <- boot(data = total, statistic = standardization, R = 10)
  
  # generating confidence intervals
  se <- c(sd(results_bootstrap$t[,1]), sd(results_bootstrap$t[,2]), 
          sd(results_bootstrap$t[,3]), sd(results_bootstrap$t[,4]))
  mean <- results_bootstrap$t0
  ll <- mean - qnorm(0.975)*se
  ul <- mean + qnorm(0.975)*se
  
  bootstrap <- data.frame(cbind(c("Observed", "No Treatment", "Treatment", 
                                  "Treatment - No Treatment"), mean, se, ll, ul)) 
  
  bootstrap <-
    bootstrap %>%
    filter(V1 == "Treatment - No Treatment") %>%
    mutate(ll = as.numeric(ll),ul = as.numeric(ul)) %>%
    mutate(coverage = if_else(ll <= ATE & ul >= ATE, 1, 0))
  
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
  ll = mean(pull(diff_results,ll)),
  ul = mean(pull(diff_results,ul)),
  coverage = mean(pull(diff_results,coverage))
)

