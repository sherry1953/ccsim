Example of a one-time simulation for a binary confounder
================
Xinyan Zhou
2021/10/25

# Binary confounder

We refer to A as the exposure, Y as the outcome, and W as the confounder
that is considered to be matched.

## Data generating mechanism

``` r
PrW <- 0.5 # Set the prevalence of confounder, W
AonY <- 0.5 # Set the strength of the exposure-disease relationship
superN <- 1000000 

WonA <- 0.388 # Set the strength of the condouner-exposure relationship
WonY <- 0.388 # Set the strength of the confounder-outcome relationship

n <- 5000 # Set the sample size

W <- rbinom(superN, 1, PrW) 
A <- rbinom(superN,1,plogis(-2.40 + WonA*W))
PrA <- mean(A) # Calculate the prevalence of the exposure
Y <- rbinom(superN,1,plogis(-2.46 + AonY*A + WonY*W))
PrY <- mean(Y) # Calculate the prevalence of the outcome
dat <- data.frame(W,A,Y) 

id <- rownames(dat)
dat <- cbind(id = id, dat)
```

## Calculate the odds ratio for the confounder-outcome association

``` r
model_WonA <- glm(formula = Y ~ W, family = binomial(),data = dat)
OR_WonA <- exp(coefficients(model_WonA))[["W"]]
```

## Calculate the true odds ratio

``` r
# define models
umodel <- "Y ~ A"
amodel <- "Y ~ A + W"

# calculate the unadjusted and adjusted ORs in the super population
unadjusted <- glm(formula = umodel,family = binomial(),data = dat)
adjusted   <- glm(formula = amodel,family = binomial(),data = dat)

uOR <- exp(coefficients(unadjusted))[["A"]]
tOR <- exp(coefficients(adjusted))[["A"]]
```

## Sampling from the super population

``` r
obsdat <- dat[sample(x = superN, size = n, replace = F), ]
```

## Matched design

In this simulation, we use frequency matching to find matched controls.

``` r
  case <- subset(obsdat, Y == 1)
  control <- subset(obsdat, Y == 0)
  
  cases_freq <- group_by(case, W) %>% 
                summarise(n = n_distinct(id))
  cs <- NA
  
  for (i in 1:nrow( cases_freq) ) {
    pts <- control[control$W == cases_freq$W[i] ,]$id 
    recs <- sample(x = 1:length(pts), size = (cases_freq$n[i])) 
    cs <- c(cs, pts[recs]) 
  }
  
  cs <- cs[2:length(cs)] 
  cs <- as.data.frame(matrix(cs, ncol = 1))
  colnames(cs)[1] <- 'id'
  
  cs <- left_join(x = cs,
                  y = control, 
                  by = ('id' = 'id'))
  
  match <- rbind(case,cs)
  
  unadjusted_m <- glm(formula = umodel,family = binomial(),data = match)
  adjusted_m   <- glm(formula = amodel,family = binomial(),data = match)
  
  uOR <- exp(coefficients(unadjusted_m))[["A"]]
  aOR <- exp(coefficients(adjusted_m))[["A"]]
  
  confuOR_lb <- exp(confint(unadjusted_m))[["A","2.5 %"]]
  confuOR_ub <- exp(confint(unadjusted_m))[["A","97.5 %"]]  
  confaOR_lb <- exp(confint(adjusted_m))[["A","2.5 %"]]
  confaOR_ub <- exp(confint(adjusted_m))[["A","97.5 %"]] 
  
  aOR_bias <- aOR - tOR
  percent_bias <- pbias(aOR,tOR)
  mean_squared_error <- mse(aOR,tOR)
  cov <- as.numeric(confaOR_lb <= tOR & confaOR_ub >= tOR)
  
  results_matched <- data.frame(cbind(uOR, aOR, 
                              uOR_lb = confuOR_lb,
                              uOR_ub = confuOR_ub,
                              aOR_lb = confaOR_lb,
                              aOR_ub = confaOR_ub,
                              range = confaOR_ub - confaOR_lb,
                              bias = aOR_bias,
                              pbias = percent_bias,
                              mse = mean_squared_error,
                              coverage = cov))
  results_matched <-
    results_matched %>%
    pivot_longer(uOR:coverage, names_to = "Measures",
               values_to = "Matched design")
```

## Unmatched design

For an unmatched design, we use random sampling to select controls.

``` r
  cases <- subset(obsdat, Y > 0.5, select = W:Y)
  controls <- subset(obsdat, Y < 0.5, select = W:Y)
  
  index <- sample(1:nrow(controls), nrow(cases),replace = FALSE)
  unmatched <- controls[index, ]
  total <- rbind(cases,unmatched)
  
  unadjusted_u <- glm(formula = umodel,family = binomial(),data = total)
  adjusted_u   <- glm(formula = amodel,family = binomial(),data = total)
  
  uOR <- exp(coefficients(unadjusted_u))[["A"]]
  aOR <- exp(coefficients(adjusted_u))[["A"]]
  
  confuOR_lb <- exp(confint(unadjusted_u))[["A","2.5 %"]]
  confuOR_ub <- exp(confint(unadjusted_u))[["A","97.5 %"]]  
  confaOR_lb <- exp(confint(adjusted_u))[["A","2.5 %"]]
  confaOR_ub <- exp(confint(adjusted_u))[["A","97.5 %"]] 
  
  aOR_bias <- aOR - tOR
  percent_bias <- pbias(aOR,tOR)
  mean_squared_error <- mse(aOR,tOR)
  cov <- as.numeric(confaOR_lb <= tOR & confaOR_ub >= tOR)
  
  results_unmatched <- data.frame(cbind(uOR, aOR, 
                              uOR_lb = confuOR_lb,
                              uOR_ub = confuOR_ub,
                              aOR_lb = confaOR_lb,
                              aOR_ub = confaOR_ub,
                              range = confaOR_ub - confaOR_lb,
                              bias = aOR_bias,
                              pbias = percent_bias,
                              mse = mean_squared_error,
                              coverage = cov))
  results_unmatched <-
    results_unmatched %>%
    pivot_longer(uOR:coverage, names_to = "Measures",
               values_to = "Unmatched design")
```

## Combine results from matched and unmatched designs

``` r
  "True odds ratio"
```

    ## [1] "True odds ratio"

``` r
  tOR
```

    ## [1] 1.672961

``` r
  left_join(results_matched, results_unmatched, by = "Measures") %>%
  knitr::kable(digits = 2)
```

| Measures | Matched design | Unmatched design |
|:---------|---------------:|-----------------:|
| uOR      |           1.63 |             2.21 |
| aOR      |           1.63 |             2.12 |
| uOR_lb   |           1.14 |             1.51 |
| uOR_ub   |           2.34 |             3.28 |
| aOR_lb   |           1.14 |             1.44 |
| aOR_ub   |           2.34 |             3.14 |
| range    |           1.20 |             1.70 |
| bias     |          -0.04 |             0.44 |
| pbias    |          -2.50 |            26.50 |
| mse      |           0.00 |             0.20 |
| coverage |           1.00 |             1.00 |
