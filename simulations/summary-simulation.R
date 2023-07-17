rm(list = ls())
library(tidyverse)
library(xtable)
summarize_sim <- function(results, Delta_C, Delta_I, m){
  results[is.infinite(results)] <- NA
  results[abs(results) > 50] <- NA
  est_index_C <- c(1,3,5,7,9,21,25)
  est_index_I <- c(11,13,15,17,19,23,27)
  sd_index_C <- c(2,4,6,8,10,22,26)
  sd_index_I <- c(12,14,16,18,20,24,28)
  scaling <- sqrt(m)/sqrt(m)
  summary_results <- cbind(deltaC_bias = apply(results[est_index_C,], 1, mean, na.rm = T) - Delta_C,
                           deltaC_emp_se = apply(results[est_index_C,], 1, sd, na.rm = T) ,
                           deltaC_avg_se = apply(results[sd_index_C,], 1, mean, na.rm = T)*scaling,
                           deltaC_cp = apply(abs((results[est_index_C,]-Delta_C)/results[sd_index_C,]/scaling) <= qt(0.975, m-5), 1, mean, na.rm = T),
                           deltaI_est = apply(results[est_index_I,], 1, mean, na.rm = T) - Delta_I,
                           deltaI_emp_se = apply(results[est_index_I,], 1, sd, na.rm = T),
                           deltaI_avg_se = apply(results[sd_index_I,], 1, mean, na.rm = T)*scaling,
                           deltaI_cp = apply(abs((results[est_index_I,]-Delta_I)/results[sd_index_I,]/scaling) <= qt(0.975, m-5), 1, mean, na.rm = T))
  rownames(summary_results) <- c("unadj", "gee", "glmm", "aug-gee", "tmle", "eff-pm", "eff-ml")
  round(summary_results,3)
}

# simulation 1 ----------------------------------
# simulation 1-1 scenario 1: m = 30 and non-informative
summarize_sim(readRDS("simulations/simulation1-1-m30-non-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=30) %>% xtable

# simulation 1-1 scenario 2: m = 30 and informative
summarize_sim(readRDS("simulations/simulation1-1-m30-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=30) %>% xtable

# simulation 1-1 scenario 3: m = 100 and non-informative
summarize_sim(readRDS("simulations/simulation1-1-m100-non-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=100) %>% xtable

# simulation 1-1 scenario 4: m = 100 and informative
summarize_sim(readRDS("simulations/simulation1-1-m100-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=100) %>% xtable


# simulation 1-2 scenario 1: m = 30 and non-informative
summarize_sim(readRDS("simulations/simulation1-2-m30-non-informative.rds"), Delta_C = 5.5, Delta_I = mean((10:100)^2)/10/55, m=30) %>% xtable

# simulation 1-2 scenario 2: m = 30 and informative
summarize_sim(readRDS("simulations/simulation1-2-m30-informative.rds"), Delta_C = 5.5, Delta_I = mean((10:100)^2)/10/55, m=30) %>% xtable

# simulation 1-2 scenario 3: m = 100 and non-informative
summarize_sim(readRDS("simulations/simulation1-2-m100-non-informative.rds"), Delta_C = 5.5, Delta_I = mean((10:100)^2)/10/55, m=100) %>% xtable

# simulation 1-2 scenario 4: m = 100 and informative
summarize_sim(readRDS("simulations/simulation1-2-m100-informative.rds"), Delta_C = 5.5, Delta_I = mean((10:100)^2)/10/55, m=100) %>% xtable


# simulation2 ----------------------------------
# simulation 2-1 scenario 1: m = 30 and non-informative
summarize_sim(readRDS("simulations/simulation2-1-m30-non-informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=30) %>% xtable

# simulation 2-1 scenario 2: m = 30 and informative
summarize_sim(readRDS("simulation2-m=30informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=30) %>% xtable

# simulation 2-1 scenario 3: m = 100 and non-informative
summarize_sim(readRDS("simulations/simulation2-1-m100-non-informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=100) %>% xtable

  # simulation 2-1 scenario 4: m = 100 and informative
summarize_sim(readRDS("simulations/simulation2-1-m100-informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=100) %>% xtable


# simulation 2-2 scenario 1: m = 30 and non-informative
summarize_sim(readRDS("simulations/simulation2-2-m30-non-informative.rds"), Delta_C = 1.55687, Delta_I = 1.348552, m=30) %>% xtable

# simulation 2-2 scenario 2: m = 30 and informative
summarize_sim(readRDS("simulations/simulation2-2-m30-informative.rds"), Delta_C = 1.55687, Delta_I = 1.348552, m=30) %>% xtable

# simulation 2-2 scenario 3: m = 100 and non-informative
summarize_sim(readRDS("simulations/simulation2-2-m100-non-informative.rds"), Delta_C = 1.55687, Delta_I = 1.348552, m=100) %>% xtable

# simulation 2-2 scenario 4: m = 100 and informative
summarize_sim(readRDS("simulations/simulation2-2-m100-informative.rds"), Delta_C = 1.55687, Delta_I = 1.348552, m=100) %>% xtable


# Additional simulations -----------------------
# simulation 1-1: m = 1000 and non-informative
summarize_sim(readRDS("simulations/simulation1-1-m1000-non-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=1000) %>% xtable

# simulation 1-1: m = 1000 and informative
summarize_sim(readRDS("simulations/simulation1-1-m1000-informative.rds"), Delta_C = 6, Delta_I = 8.666667, m=1000) %>% xtable

# simulation 2-1: m = 1000 and non-informative
summarize_sim(readRDS("simulations/simulation2-1-m1000-non-informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=1000) %>% xtable

# simulation 2-1: m = 1000 and informative
summarize_sim(readRDS("simulations/simulation2-1-m1000-informative.rds"), Delta_C = 1.538776, Delta_I = 1.184246, m=1000) %>% xtable


