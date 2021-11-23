library(hydroGOF)
library(tidyverse)
library(splitstackshape)

rm(list = ls())

PrW <- 0.5
AonY <- 0.5
superN <- 1000000
simN <- 5

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
  
  cs <- left_join(x = cs,
                  y = controls, 
                  by = ('id' = 'id'))
  
  match <- rbind(cases,cs)
  
  ## Create a pseudo cohort from the case-control sample
  
  cs_W0_A <- 
    cs %>%
    filter(W == 0)
  
  prop_W0_A <-
    cs_W0_A %>%
    group_by(A) %>%
    summarise(prop = n()/nrow(cs_W0_A))
  
  prop_W0_A <- pull(prop_W0_A, prop)
  prop_W0_A0 <- prop_W0_A[1]
  prop_W0_A1 <- prop_W0_A[2]
  
  cs_W1_A <- 
    cs %>%
    filter(W == 1)
  
  prop_W1_A <-
    cs_W1_A %>%
    group_by(A) %>%
    summarise(prop = n()/nrow(cs_W1_A))
  
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

diff_list <- lapply(1:simN, function(x) matched(obsdat))
diff_results <- do.call(rbind, diff_list)

final_results <- tibble(
  ATE,
  RD_cohort = mean(pull(diff_results,RD_cohort)),
  RD_step1 = mean(pull(diff_results,RD_step1)),
  RD_step2 = mean(pull(diff_results,RD_step2)),
  diff_1 = mean(pull(diff_results,diff_1)),
  diff_2 = mean(pull(diff_results,diff_2))
)
