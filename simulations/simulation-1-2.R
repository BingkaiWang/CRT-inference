# ------------
# simulation 1 (continuous outcomes)
# Estimand: average treatment effect with Delta_C = 5.5 and Delta_I = 6.75
# Methods: unadjusted, GEE, LMM, Aug-GEE, TMLE, Eff-PM, Eff-ML
# The population size for each cluster follows Unif[10,100], which is different from simulation-1-1.
# ------------
rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack) #GEE
library(lme4) # GLMM
library(CRTgeeDR) #Aug-GEE
library(SuperLearner) # TMLE
library(tmle) # TMLE
library(foreach)
library(doSNOW)
cl <- makeCluster(16)
registerDoSNOW(cl)
expit <- function(x){1/(1+exp(-x))}

# simulation ----------------
m <- 30 # number of clusters
pi <- 0.5 # randomization ratio
sim_size <- 10000
package_list <- c("tidyverse", "lme4", "geepack", "CRTgeeDR", "tmle", "SuperLearner")
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
simulation_scenario <- "informative" # or "non-informative"
print(paste0("simulation1-2-m",m,"-",simulation_scenario))

tictoc::tic()
results <- foreach(iter = 1:sim_size, .combine = cbind, .packages = package_list) %dopar% {
  N <- sample(10:100, size = m, replace = T)
  C1 <- rnorm(m, mean = 0.1 * N/2, sd = 2)
  C2 <- rbinom(m, size = 1, prob = expit(C1 * log(N/20)))
  A <- rbinom(m, size = 1, prob = pi)
  complete_data <- map(1:m, function(j) {
    X1 <- rbinom(N[j], size = 1, prob = N[j]/100)
    X2 <- rnorm(N[j], mean = ifelse(C2[j], 1, -1) * mean(X1), sd = 3)
    if(simulation_scenario == "non-informative"){
      M1 <- 10 + sample(c(-1,0), size = 1)
      M0 <- 10 + sample(c(-1,0), size = 1)
    } else if(simulation_scenario == "informative"){
      M1 <- round(N[j]/10) + 5 * C2[j]
      M0 <- ifelse(N[j] > 55, 6, 3)
    }
    Y1 <- sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 + 5 * exp(X1) * abs(X2) + rnorm(N[j], sd = 1) + N[j]/10
    Y0 <- sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 + 5 * exp(X1) * abs(X2) + rnorm(N[j], sd = 1) + rnorm(1, mean = 0, sd = 1)
    cluster_id <- j
    data.frame(cbind(Y1, Y0, M1, M0, X1, X2, A=A[j], N=N[j], C1=C1[j], C2=C2[j]), cluster_id)
  })
  
  observed_data <- map_dfr(complete_data, function(d){
    AA <- d$A[1]
    MM <- d$M1[1] * AA + d$M0[1] * (1-AA)
    NN <- d$N[1]
    obs_indi <- sample(1:NN, size = MM, replace = F)
    d[obs_indi,] %>% 
      mutate(Y = AA * Y1 + (1-AA) * Y0, M = MM) %>%
      dplyr::select(-(Y1:M0))
  })
  cl_data <- observed_data %>% group_by(cluster_id) %>% summarise_all(mean)
  
  tryCatch({
    # estimating Delta_C ------
    unadj_est <- mean(cl_data$Y[cl_data$A==1]) - mean(cl_data$Y[cl_data$A==0])
    unadj_se <- sqrt(var(cl_data$Y[cl_data$A==1])/pi + var(cl_data$Y[cl_data$A==0])/(1-pi))/sqrt(m-1)
    
    gee <- geeglm(Y~ A + N + X1 + X2 + C1 + C2, id = cluster_id, 
                  data =observed_data, family = gaussian, corstr = "exchangeable")
    gee_est <- summary(gee)$coefficients['A', 'Estimate']
    gee_se <- summary(gee)$coefficients['A', 'Std.err']
    
    glmm <- lmer(Y~(1|cluster_id) + A + N + X1 + X2 + C1 + C2, data = observed_data)
    glmm_est <- summary(glmm)$coefficients['A', 'Estimate']
    glmm_se <- summary(glmm)$coefficients['A', 'Std. Error']
    
    aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = observed_data, family = "gaussian",corstr = "exchangeable",
                               nameY = "Y", nameTRT = "A", model.augmentation.trt = Y ~ N + X1 + X2 + C1 + C2,
                               model.augmentation.ctrl = Y ~ N + X1 + X2 + C1 + C2)
    aug_gee_est <- summary(aug_gee)$beta[2]
    aug_gee_se <- summary(aug_gee)$se.robust[2]
    
    cl_tmle <- tmle(Y = cl_data$Y, A = cl_data$A, W = cl_data[,c('N', 'C1', 'C2')], 
                    Q.SL.library = SLmethods, g1W = rep(pi, m), family = "gaussian")
    cl_est <- cl_tmle$estimates$ATE$psi
    cl_se <- sqrt(cl_tmle$estimates$ATE$var.psi)
    
    result_C <- c(unadj_est, unadj_se, gee_est, gee_se, glmm_est, glmm_se, aug_gee_est, aug_gee_se, cl_est, cl_se)
    
    # estimating Delta_I ------
    unadj_mu1 <- sum(cl_data$Y * cl_data$A * cl_data$N)/sum(cl_data$A * cl_data$N)
    unadj_mu0 <- sum(cl_data$Y * (1-cl_data$A) * cl_data$N)/sum((1-cl_data$A) * cl_data$N)
    unadj_est <- unadj_mu1 - unadj_mu0
    unadj_se <- sd(cl_data$N/mean(cl_data$N) * (cl_data$A*(cl_data$Y-unadj_mu1)/pi - (1-cl_data$A)*(cl_data$Y - unadj_mu0)/(1-pi)))/sqrt(m-1)
    
    gee <- geeglm(Y~ A + N + X1 + X2 + C1 + C2, id = cluster_id, weights = observed_data$N,
                  data =observed_data, family = gaussian, corstr = "exchangeable")
    gee_est <- summary(gee)$coefficients['A', 'Estimate']
    gee_se <- summary(gee)$coefficients['A', 'Std.err']
    
    scaled_data <- map_dfr(1:m, function(j){
      dd <- observed_data[observed_data$cluster_id == j,]
      N_j <- dd$N[1]
      map_dfr(1: N_j, function(k){
        mutate(dd, cluster_id = as.numeric(paste0(j, k)))
      })
    })
    glmm <- lmer(Y~(1|cluster_id) + A + N + X1 + X2 + C1 + C2, REML=F, data = scaled_data)
    glmm_est <- summary(glmm)$coefficients['A', 'Estimate']
    glmm_se <- summary(glmm)$coefficients['A', 'Std. Error'] * sqrt(nrow(scaled_data)) / sqrt(nrow(observed_data))
    
    aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = scaled_data, family = "gaussian",corstr = "exchangeable",
                               nameY = "Y", nameTRT = "A", model.augmentation.trt = Y ~ N + X1 + X2 + C1 + C2,
                               model.augmentation.ctrl = Y ~ N + X1 + X2 + C1 + C2)
    aug_gee_est <- summary(aug_gee)$beta[2]
    aug_gee_se <- summary(aug_gee)$se.robust[2] * sqrt(nrow(scaled_data)) / sqrt(nrow(observed_data))
    
    I_tmle <- tmle(Y = observed_data$Y, A = observed_data$A, W = observed_data[,c('N', 'X1', 'X2', 'C1', 'C2')], 
                   Q.SL.library = SLmethods, g1W = rep(pi, nrow(observed_data)), family = "gaussian")
    tmle_est <- I_tmle$estimates$ATE$psi
    tmle_se <- cbind(IC = I_tmle$estimates$IC$IC.ATE, N = observed_data$N, id = observed_data$cluster_id) %>%
      as.data.frame %>% group_by(id) %>% summarise_all(sum) %>% dplyr::select(IC) %>% var %>% sqrt(.)/ nrow(observed_data) * sqrt(m-5)
    
    
    result_I <- c(unadj_est, unadj_se, gee_est, gee_se, glmm_est, glmm_se, aug_gee_est, aug_gee_se, tmle_est, tmle_se)
    
    # proposed method -----
    dd <- left_join(observed_data, cl_data[,c("cluster_id", "X1", "X2")], by = "cluster_id")
    # eff-pm -----
    kappa1 <- glm(A~M+N+C1+C2, family = binomial, data = cl_data) %>% predict(type = "response")
    kappa0 <- 1 - kappa1
    zeta1 <- predict(lm(Y~N+C1+C2, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
    zeta0 <- predict(lm(Y~N+C1+C2, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
    eta <- dd %>% mutate(eta1= predict(lm(Y~., data = dplyr::select(dd[dd$A==1,],-c('A', 'cluster_id'))), newdata = dd, type = "response")) %>%
      mutate(eta0 = predict(lm(Y~., data = dplyr::select(dd[dd$A==0,],-c('A', 'cluster_id'))), newdata = mutate(dd, A = 0), type = "response")) %>%
      dplyr::select(cluster_id, eta1, eta0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean) %>%
      as.data.frame()
    mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
    mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
    eff_pm_est_C <- mean(mu_C1) - mean(mu_C0)
    eff_pm_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
    eff_pm_se_C <- sqrt(var(mu_C1 - mu_C0)/(m-5))
    eff_pm_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - eff_pm_est_I)/mean(cl_data$N))/(m-5))
    
    # eff-ML -----
    kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("M", "N", "C1", "C2")], family = "binomial",
                           SL.library = SLmethods, cvControl = list(V=5))$SL.predict
    kappa0 <- 1 - kappa1
    zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N", "C1", "C2")], family = "gaussian", 
                              SL.library = SLmethods, cvControl = list(V=5))
    zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N", "C1", "C2")], family = "gaussian", 
                              SL.library = SLmethods, cvControl = list(V=5))
    zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N", "C1", "C2")])$pred
    zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N", "C1", "C2")])$pred
    eta_fit1 <- SuperLearner(dd$Y[dd$A==1], X = dplyr::select(dd[dd$A==1,], -c("Y", "cluster_id","A")), family = "gaussian", 
                             id = dd$cluster_id[dd$A==1], SL.library = SLmethods, cvControl = list(V=5))
    eta_fit0 <- SuperLearner(dd$Y[dd$A==0], X = dplyr::select(dd[dd$A==0,], -c("Y", "cluster_id","A")), family = "gaussian", 
                             id = dd$cluster_id[dd$A==0], SL.library = SLmethods, cvControl = list(V=5))
    eta <- dd %>% mutate(eta1= predict(eta_fit1, newdata = dplyr::select(dd, -c("Y", "cluster_id","A")) %>% mutate(A = 1))$pred) %>%
      mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(dd, -c("Y", "cluster_id",'A')) %>% mutate(A = 0))$pred) %>%
      dplyr::select(cluster_id, eta1, eta0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean) %>%
      as.data.frame()
    mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
    mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
    eff_ml_est_C <- mean(mu_C1) - mean(mu_C0)
    eff_ml_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
    eff_ml_se_C <- sqrt(var(mu_C1 - mu_C0)/(m-5))
    eff_ml_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - eff_pm_est_I)/mean(cl_data$N))/(m-5))
    result_eff <- c(eff_pm_est_C, eff_pm_se_C, eff_pm_est_I, eff_pm_se_I, eff_ml_est_C, eff_ml_se_C, eff_ml_est_I, eff_ml_se_I)
    
    c(result_C, result_I, result_eff)
  }, error = function(e) { return(rep(NA, 28))})
  
  
}
tictoc::toc()

delta_C <- 5.5
delta_I <- mean((10:100)^2)/10/55
stopCluster(cl)

est_index_C <- c(1,3,5,7,9,21,25)
est_index_I <- c(11,13,15,17,19,23,27)
sd_index_C <- c(2,4,6,8,10,22,26)
sd_index_I <- c(12,14,16,18,20,24,28)
summary_results <- cbind(deltaC_bias = apply(results[est_index_C,], 1, mean, na.rm = T) - delta_C,
                         deltaC_emp_se = apply(results[est_index_C,], 1, sd, na.rm = T),
                         deltaC_avg_se = apply(results[sd_index_C,], 1, mean, na.rm = T),
                         deltaC_cp = apply(abs((results[est_index_C,]-delta_C)/results[sd_index_C,]) <= qt(0.975, m), 1, mean, na.rm = T),
                         deltaI_bias = apply(results[est_index_I,], 1, mean, na.rm = T) - delta_I,
                         deltaI_emp_se = apply(results[est_index_I,], 1, sd, na.rm = T),
                         deltaI_avg_se = apply(results[sd_index_I,], 1, mean, na.rm = T),
                         deltaI_cp = apply(abs((results[est_index_I,]-delta_I)/results[sd_index_I,]) <= qt(0.975, m), 1, mean, na.rm = T))

rownames(summary_results) <- c("unadj", "gee", "glmm", "aug-gee", "tmle", "eff-pm", "eff-ml")
saveRDS(results, paste0("simulation1-2-m",m,"-",simulation_scenario,".rds"))
round(summary_results,2)
