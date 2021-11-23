library(hydroGOF)
library(tidyverse)
library(splitstackshape)

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
Y <- rbinom(superN,1,plogis(-3.46 + AonY*A + WonY*W))
PrY <- mean(Y)
Y1 <- rbinom(superN,1,plogis(-2.46 + AonY + WonY*W))
Y0 <- rbinom(superN,1,plogis(-2.46 + WonY*W))
R1 <- mean(Y1)
R0 <- mean(Y0)
ATE <- R1 - R0

dat <- data.frame(W,A,Y)
id <- rownames(dat)
dat <- cbind(id = id, dat)

unmatched <- function(obsdat){
  
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
  
  cases <- subset(obsdat, Y == 1, select = W:Y)
  controls <- subset(obsdat, Y == 0, select = W:Y)
  
  prop_W <- 
    controls %>%
    group_by(W) %>%
    summarise(prop = n()/nrow(controls))
  
  prop_W <- pull(prop_W, prop)
  prop_W0 <- prop_W[1]
  prop_W1 <- prop_W[2]
  
  ## Select controls
  
  index <- sample(1:nrow(controls), nrow(cases),replace = FALSE)
  controls_selected <- controls[index, ]
  
  ## Create a pseudo cohort from the case-control sample
  
  type_prop_step1 <-
    controls_selected %>%
    mutate(type = if_else(W == 0 & A == 0, 1, 
                          if_else(W == 0 & A == 1, 2,
                                  if_else(W == 1 & A == 0, 3,
                                          if_else(W == 1 & A == 1, 4, 0))))) %>%
    count(type) %>%
    mutate(prop = n/nrow(controls_selected))
  
  prop_step1 <- pull(type_prop_step1,prop)
  prop_type_1 <- prop_step1[1]
  prop_type_2 <- prop_step1[2]
  prop_type_3 <- prop_step1[3]
  prop_type_4 <- prop_step1[4]
  
  type_1 <-
    tibble(W = 0, A = 0, Y = 0) %>%
    expandRows(type_1, count = prop_type_1 * nrow(controls) , count.is.col = FALSE)
  
  type_2 <-
    tibble(W = 0, A = 1, Y = 0) %>%
    expandRows(type_2, count = prop_type_2 * nrow(controls) , count.is.col = FALSE)
  
  type_3 <-
    tibble(W = 1, A = 0, Y = 0) %>%
    expandRows(type_3, count = prop_type_3 * nrow(controls) , count.is.col = FALSE)
  
  type_4 <-
    tibble(W = 1, A = 1, Y = 0) %>%
    expandRows(type_4, count = prop_type_4 * nrow(controls) , count.is.col = FALSE)
  
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
  
  ## Standardization
  
  total_type <-
    total %>%
    mutate(type = if_else(W == 0 & A == 0, 1, 
                  if_else(W == 0 & A == 1, 2,
                  if_else(W == 1 & A == 0, 3,
                  if_else(W == 1 & A == 1, 4, 0))))) 
  
  type_1_step2 <-
    total_type %>%
    filter(type == 1)
  E_Y1_type_1 <- mean(pull(type_1_step2,Y))
  
  type_2_step2 <-
    total_type %>%
    filter(type == 2)
  E_Y1_type_2 <- mean(pull(type_2_step2,Y))
  
  type_3_step2 <-
    total_type %>%
    filter(type == 3)
  E_Y1_type_3 <- mean(pull(type_3_step2,Y))
  
  type_4_step2 <-
    total_type %>%
    filter(type == 4)
  E_Y1_type_4 <- mean(pull(type_4_step2,Y))
  
  df_PrW <-
    total %>%
    count(W) %>%
    mutate(prop = n/nrow(total))
  PrW0 <- df_PrW[1,3]
  PrW1 <- df_PrW[2,3]
  
  R1_step2 <-  E_Y1_type_2*PrW0 + E_Y1_type_4*PrW1
  R0_step2 <-  E_Y1_type_1*PrW0 + E_Y1_type_3*PrW1 
  
  RD_step2 <- R1_step2 - R0_step2
  
  diff_1 <- RD_cohort - RD_step1
  diff_2 <- ATE - RD_step2
  
  results <-
    tibble(RD_cohort,RD_step1,RD_step2,diff_1,diff_2)
  
  return(results)
}

diff_list <- lapply(1:simN, function(x) unmatched(obsdat))
diff_results <- do.call(rbind, diff_list)

final_results <- tibble(
  ATE,
  RD_cohort = mean(pull(diff_results,RD_cohort)),
  RD_step1 = mean(pull(diff_results,RD_step1)),
  RD_step2 = mean(pull(diff_results,RD_step2)),
  diff_1 = mean(pull(diff_results,diff_1)),
  diff_2 = mean(pull(diff_results,diff_2))
)
