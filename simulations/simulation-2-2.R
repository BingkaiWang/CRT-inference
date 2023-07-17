# ------------
# simulation 2 (binary outcomes)
# Estimand: risk ratio with Delta_C = 0.79 and Delta_I = 0.71
# Methods to compare: GEE, LMM, Aug-GEE, cluster-level ANCOVA, individual-level ANCOVA
# ------------
rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack) #GEE
library(lme4) # GLMM
library(CRTgeeDR) #Aug-GEE
# library(ltmle) # TMLE
library(SuperLearner) # TMLE
library(tmle) # TMLE
library(foreach)
library(doSNOW)
cl <- makeCluster(16)
registerDoSNOW(cl)
# source("Stage2_Function.R") # for TMLE
expit <- function(x){1/(1+exp(-x))}
var_rr <- function(IF, mu){
  c(1/mu[2], -mu[1]/mu[2]^2) %*% var(IF) %*% c(1/mu[2], -mu[1]/mu[2]^2)
}
var_aug_gee <- function(beta, var){
  derivative <- c(-exp(-beta[1]) * (1-exp(-beta[2])), (1+exp(-beta[1])) * exp(-beta[1]-beta[2])) * expit(beta[1] + beta[2])^2
  derivative %*% var %*% derivative
}
simulation_scenario <- "non-informative" # or "non-informative"

# find true Delta_C and Delta_I -------------
m<- 10000
pi <- 0.5
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
  Y1 <- rbinom(N[j], size = 1, prob = expit(sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 +   1.5 * exp(X1) * sqrt(abs(X2)) - N[j]/40))
  Y0 <- rbinom(N[j], size = 1, prob = expit(sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 +   1.5 * (2 *X1 - 1)  * sqrt(abs(X2)) + rnorm(N[j], sd = 1)))
  cluster_id <- j
  data.frame(cbind(Y1, Y0, M1, M0, X1, X2, A=A[j], N=N[j], C1=C1[j], C2=C2[j]), cluster_id)
})  
delta_C <- map_dfr(complete_data, ~c(y1 = mean(.$Y1), y0 = mean(.$Y0))) %>% apply(2, mean) %>% {.[1]/.[2]}
delta_I <- map_dfr(complete_data, ~c(y1 = sum(.$Y1), y0 = sum(.$Y0))) %>% apply(2, mean) %>% {.[1]/.[2]}
delta_C
delta_I

# these are delta_C and delta_I for the OBSERVED people (instead of the source population), which can be viewed as the estimand for GEE
(map_dfr(complete_data, ~c(y1 = mean(.$Y1) *.$M1[1], y0 = mean(.$Y0) * .$M0[1])) %>% apply(2, mean) %>% {.[1]/.[2]}) /
  (map_dfr(complete_data, ~c(y1 = .$M1[1], y0 =  .$M0[1])) %>% apply(2, mean) %>% {.[1]/.[2]})
(map_dfr(complete_data, ~c(y1 = sum(.$Y1) *.$M1[1], y0 = sum(.$Y0) * .$M0[1])) %>% apply(2, mean) %>% {.[1]/.[2]}) /
  (map_dfr(complete_data, ~c(y1 = .$M1[1] * .$N[1], y0 =  .$M0[1] * .$N[1])) %>% apply(2, mean) %>% {.[1]/.[2]})



# simulation ------------
m <- 30 # number of clusters
pi <- 0.5 # randomization ratio
sim_size <- 10000
package_list <- c("tidyverse", "lme4", "geepack", "CRTgeeDR", "tmle", "SuperLearner")
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
print(paste0("simulation2-2-m",m,"-",simulation_scenario))

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
    Y1 <- rbinom(N[j], size = 1, prob = expit(sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 +   1.5 * exp(X1) * sqrt(abs(X2)) - N[j]/40))
    Y0 <- rbinom(N[j], size = 1, prob = expit(sin(C1[j]) * (2 *C2[j] - 1) * N[j]/60 +   1.5 * (2 *X1 - 1)  * sqrt(abs(X2)) + rnorm(N[j], sd = 1)))
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
    # Unadjusted -----
    mu_1 <- mean(cl_data$Y[cl_data$A==1])
    mu_0 <- mean(cl_data$Y[cl_data$A==0])
    unadj_est <- mu_1/mu_0
    unadj_var <- var(cl_data$Y[cl_data$A==1])/pi/mu_0^2 + var(cl_data$Y[cl_data$A==0])/(1-pi)/mu_0^4*mu_1^2
    unadj_se <- sqrt(unadj_var)/sqrt(m-1)
    # gee ----
    gee <- geeglm(Y~ A + N + X1 + X2 + C1 + C2, id = cluster_id, 
                  data =observed_data, family = 'binomial', corstr = "independence")
    gee_data <- observed_data %>% dplyr::select(cluster_id, A, Y) %>%
      mutate(pred1 = predict(gee, newdata = mutate(observed_data, A = 1), type = 'response')) %>%
      mutate(pred0 = predict(gee, newdata = mutate(observed_data, A = 0), type = 'response')) %>%
      group_by(cluster_id) %>% summarise(M = n(), A = mean(A), Y = mean(Y), pred1 = mean(pred1), pred0 = mean(pred0))
    gee_muC_hat <- colMeans(gee_data[,c('pred1', 'pred0')])
    gee_IF <- gee_data %>% 
      mutate(IF1 = A/pi/mean(M[gee_data$A==1]) * M * (Y-pred1) + pred1 - mean(pred1)) %>%
      mutate(IF0 = (1-A)/(1-pi)/mean(M[gee_data$A==0]) * M * (Y-pred0) + pred0 - mean(pred0)) %>%
      dplyr::select(IF1, IF0)
    gee_est <- gee_muC_hat[1]/gee_muC_hat[2]
    gee_se <- sqrt(var_rr(gee_IF,gee_muC_hat)/(m-5))
    # glmm -----
    glmm <- lmer(Y~(1|cluster_id) + A + N + X1 + X2 + C1 + C2, data = observed_data)
    glmm_data <- observed_data %>% dplyr::select(cluster_id, A, Y) %>%
      mutate(pred1 = predict(glmm, newdata = mutate(observed_data, A = 1), re.form = NA, type = 'response')) %>%
      mutate(pred0 = predict(glmm, newdata = mutate(observed_data, A = 0), re.form = NA,type = 'response')) %>%
      group_by(cluster_id) %>% summarise(M = n(), A = mean(A), Y = mean(Y), pred1 = mean(pred1), pred0 = mean(pred0)) 
    tauhat <- summary(glmm)$varcor$cluster_id[1]
    sigmahat <- summary(glmm)$sigma^2
    glmm_IF <- glmm_data %>%
      mutate(IF1 = A/pi/mean(M[glmm_data$A==1]/(sigmahat + tauhat * M[glmm_data$A==1]))/(sigmahat + tauhat * M) * M * (Y-pred1) + pred1 - mean(pred1)) %>%
      mutate(IF0 = (1-A)/(1-pi)/mean(M[glmm_data$A==0]/(sigmahat + tauhat * M[glmm_data$A==0]))/(sigmahat + tauhat * M) * M * (Y-pred0) + pred0 - mean(pred0)) %>%
      dplyr::select(IF1, IF0)
    glmm_muC_hat <- colMeans(glmm_data[,c('pred1', 'pred0')])
    glmm_est <-  glmm_muC_hat[1]/ glmm_muC_hat[2]
    glmm_se <- sqrt(var_rr(glmm_IF,glmm_muC_hat)/(m-5))
    # aug-gee -----
    aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = observed_data, family = "binomial",corstr = "exchangeable",
                               nameY = "Y", nameTRT = "A", model.augmentation.trt = Y ~ N + X1 + X2 + C1 + C2,
                               model.augmentation.ctrl = Y ~ N + X1 + X2 + C1 + C2)
    aug_gee_est <- expit(aug_gee$beta[1] + aug_gee$beta[2])/expit(aug_gee$beta[1])
    aug_gee_se <- sqrt(var_aug_gee(aug_gee$beta, aug_gee$var))
    
    # TMLE -----
    cl_tmle <- tmle(Y = cl_data$Y, A = cl_data$A, W = cl_data[,c('N', 'C1', 'C2')],
                    Q.SL.library = SLmethods, g1W = rep(pi, m), family = "gaussian")
    tmle_pred1 <- cl_tmle$Qstar[,2]
    tmle_pred0 <- cl_tmle$Qstar[,1]
    cl_est <- mean(tmle_pred1)/mean(tmle_pred0)
    tmle_IF1 <- cl_data$A / pi * (cl_data$Y - tmle_pred1) + tmle_pred1 - mean(tmle_pred1)
    tmle_IF0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - tmle_pred0) + tmle_pred0 - mean(tmle_pred0)
    cl_se <- sqrt(var_rr(cbind(tmle_IF1, tmle_IF0), c(mean(tmle_pred1), mean(tmle_pred0)))/(m-5))
    
    # results_C ---------
    result_C <- c(unadj_est, unadj_se, gee_est, gee_se, glmm_est, glmm_se, aug_gee_est, aug_gee_se, cl_est, cl_se)
    
    # estimating Delta_I ------
    scaled_data <- map_dfr(1:m, function(j){
      dd <- observed_data[observed_data$cluster_id == j,]
      N_j <- dd$N[1]
      map_dfr(1: N_j, function(k){
        mutate(dd, cluster_id = as.numeric(paste0(j, k)))
      })
    })
    # Unadjusted -----
    unadj_mu1 <- sum(cl_data$Y * cl_data$A * cl_data$N)/sum(cl_data$A * cl_data$N)
    unadj_mu0 <- sum(cl_data$Y * (1-cl_data$A) * cl_data$N)/sum((1-cl_data$A) * cl_data$N)
    unadj_est <- unadj_mu1/unadj_mu0
    unadj_IF1 <- cl_data$N/mean(cl_data$N) * cl_data$A/pi * (cl_data$Y-unadj_mu1)
    unadj_IF0 <- cl_data$N/mean(cl_data$N) * (1-cl_data$A)/(1-pi) * (cl_data$Y - unadj_mu0)
    unadj_se <- sqrt(var_rr(cbind(unadj_IF1, unadj_IF0),c(unadj_mu1, unadj_mu0))/(m-1))
    # gee ------
    gee <- geeglm(Y~ A + N + X1 + X2 + C1 + C2, id = cluster_id,
                  data =scaled_data, family = 'binomial', corstr = "independence")
    gee_data <- observed_data %>% dplyr::select(cluster_id, A, Y, N) %>%
      mutate(pred1 = predict(gee, newdata = mutate(observed_data, A = 1), type = 'response')) %>%
      mutate(pred0 = predict(gee, newdata = mutate(observed_data, A = 0), type = 'response')) %>%
      group_by(cluster_id) %>%
      summarise(M = n(), A = mean(A), N = mean(N), Y = mean(Y), pred1 = mean(pred1), pred0 = mean(pred0))
    gee_muI1_hat <- sum(gee_data$pred1 * gee_data$N)/sum(gee_data$N)
    gee_muI0_hat <- sum(gee_data$pred0 * gee_data$N)/sum(gee_data$N)
    gee_IF <- gee_data %>%
      mutate(IF1 = N * (A/pi/mean(M[gee_data$A==1]) * M * (Y-pred1) + pred1 - gee_muI1_hat)/mean(N)) %>%
      mutate(IF0 = N * ((1-A)/(1-pi)/mean(M[gee_data$A==0]) * M * (Y-pred0) + pred0 - gee_muI0_hat)/mean(N)) %>%
      dplyr::select(IF1, IF0)
    gee_est <- gee_muI1_hat/gee_muI0_hat
    gee_se <- sqrt(var_rr(gee_IF,c(gee_muI1_hat, gee_muI0_hat))/(m-5))
    # glmm -------
    glmm <- lmer(Y~(1|cluster_id) + A + N + X1 + X2 + C1 + C2, data = scaled_data)
    glmm_data <- observed_data %>% dplyr::select(cluster_id, A, Y, N) %>%
      mutate(pred1 = predict(glmm, newdata = mutate(observed_data, A = 1), re.form = NA, type = 'response')) %>%
      mutate(pred0 = predict(glmm, newdata = mutate(observed_data, A = 0), re.form = NA, type = 'response')) %>%
      group_by(cluster_id) %>%
      summarise(M = n(), A = mean(A), N = mean(N), Y = mean(Y), pred1 = mean(pred1), pred0 = mean(pred0)) 
    tauhat <- summary(glmm)$varcor$cluster_id[1]
    sigmahat <- summary(glmm)$sigma^2
    glmm_muI1_hat <- sum(glmm_data$pred1 * glmm_data$N)/sum(glmm_data$N)
    glmm_muI0_hat <- sum(glmm_data$pred0 * glmm_data$N)/sum(glmm_data$N)
    glmm_IF <- glmm_data %>%
      mutate(IF1 = N * (A/pi/mean(M[glmm_data$A==1]/(sigmahat + tauhat * M[glmm_data$A==1])) / (sigmahat + tauhat * M) * M * (Y-pred1) + pred1 - mean(pred1)) /mean(N))%>%
      mutate(IF0 = N * ((1-A)/(1-pi)/mean(M[glmm_data$A==0]/(sigmahat + tauhat * M[glmm_data$A==0])) / (sigmahat + tauhat * M) * M * (Y-pred0) + pred0 - mean(pred0))/mean(N)) %>%
      dplyr::select(IF1, IF0)
    glmm_est <- glmm_muI1_hat/glmm_muI0_hat
    glmm_se <- sqrt(var_rr(glmm_IF,c(glmm_muI1_hat, glmm_muI0_hat))/(m-5))
    # aug-gee -------
    aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = scaled_data, family = "binomial",corstr = "exchangeable",
                               nameY = "Y", nameTRT = "A", model.augmentation.trt = Y ~ N + X1 + X2 + C1 + C2,
                               model.augmentation.ctrl = Y ~ N + X1 + X2 + C1 + C2)
    aug_gee_est <- expit(aug_gee$beta[1] + aug_gee$beta[2])/expit(aug_gee$beta[1])
    aug_gee_se <- sqrt(var_aug_gee(aug_gee$beta, aug_gee$var) * nrow(scaled_data) / nrow(observed_data))
    # TMLE ---------
    I_tmle <- tmle(Y = observed_data$Y, A = observed_data$A, W = observed_data[,c('N', 'X1', 'X2', 'C1', 'C2')],
                   Q.SL.library = SLmethods, g1W = rep(pi, nrow(observed_data)), family = "binomial")
    tmle_est <- I_tmle$estimates$RR$psi
    tmle_se <- data.frame(id = observed_data$cluster_id, IC = I_tmle$estimates$IC$IC.logRR * exp(I_tmle$estimates$RR$log.psi)) %>%
      group_by(id) %>% summarise_all(sum) %>% dplyr::select(IC) %>% var %>% sqrt(.)/ nrow(observed_data) * sqrt(m-5)
    
    # result_I -------
    result_I <- c(unadj_est, unadj_se, gee_est, gee_se, glmm_est, glmm_se, aug_gee_est, aug_gee_se, tmle_est, tmle_se)
    
    # proposed method -----
    dd <- left_join(observed_data, cl_data[,c("cluster_id", "X1", "X2")], by = "cluster_id")
    # eff-pm -----
    kappa1 <- glm(A~M+N+C1+C2, family = binomial, data = cl_data) %>% predict(type = "response")
    # kappa1 <- pi
    # kappa1 <- abs(ifelse(cl_data$N > 49, 10, 2) + 5 * cl_data$C2 - cl_data$M) < 1e-3
    kappa0 <- 1 - kappa1
    zeta1 <- predict(lm(Y~N+C1+C2, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
    zeta0 <- predict(lm(Y~N+C1+C2, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
    eta <- dd %>% mutate(eta1= predict(glm(Y~., data = dplyr::select(dd[dd$A==1,],-c('A', 'cluster_id')), family = binomial), newdata = dd, type = "response")) %>%
      mutate(eta0 = predict(glm(Y~., data = dplyr::select(dd[dd$A==0,],-c('A', 'cluster_id')), family = binomial), newdata = dd, type = "response")) %>%
      dplyr::select(cluster_id, eta1, eta0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean) %>%
      as.data.frame()
    mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
    mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
    eff_pm_est_C <- mean(mu_C1)/mean(mu_C0)
    eff_pm_est_I <- mean(cl_data$N * mu_C1)/mean(cl_data$N * mu_C0)
    eff_pm_se_C <- sqrt(var_rr(cbind(mu_C1 - mean(mu_C1),mu_C0 - mean(mu_C0)), c(mean(mu_C1), mean(mu_C0)))/(m-5))
    eff_pm_se_I <- sqrt(var_rr(cbind(cl_data$N * (mu_C1 - mean(mu_C1))/mean(cl_data$N), cl_data$N * (mu_C0 - mean(mu_C0))/mean(cl_data$N)), 
                                 c(mean(cl_data$N * mu_C1)/mean(cl_data$N), mean(cl_data$N * mu_C0)/mean(cl_data$N)))/(m-5)) 
    
    # eff-ML -----
    kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("M", "N", "C1", "C2")], family = "binomial",
                           SL.library = SLmethods, cvControl = list(V=5))$SL.predict
    # kappa1 <- abs(ifelse(cl_data$N > 49, 10, 2) + 5 * cl_data$C2 - cl_data$M) < 1e-3
    # kappa1 <- pi
    kappa0 <- 1 - kappa1
    zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N", "C1", "C2")], family = "gaussian", 
                              SL.library = SLmethods, cvControl = list(V=5))
    zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N", "C1", "C2")], family = "gaussian", 
                              SL.library = SLmethods, cvControl = list(V=5))
    zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N", "C1", "C2")])$pred
    zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N", "C1", "C2")])$pred
    eta_fit1 <- SuperLearner(dd$Y[dd$A==1], X = dplyr::select(dd[dd$A==1,], -c("Y", "cluster_id","A")), family = "binomial", 
                             id = dd$cluster_id[dd$A==1], SL.library = SLmethods, cvControl = list(V=5))
    eta_fit0 <- SuperLearner(dd$Y[dd$A==0], X = dplyr::select(dd[dd$A==0,], -c("Y", "cluster_id","A")), family = "binomial", 
                             id = dd$cluster_id[dd$A==0], SL.library = SLmethods, cvControl = list(V=5))
    eta <- dd %>% mutate(eta1= predict(eta_fit1, newdata = dplyr::select(dd, -c("Y", "cluster_id","A")))$pred) %>%
      mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(dd, -c("Y", "cluster_id",'A')))$pred) %>%
      dplyr::select(cluster_id, eta1, eta0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean) %>%
      as.data.frame()
    mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
    mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
    eff_ml_est_C <- mean(mu_C1)/mean(mu_C0)
    eff_ml_est_I <- mean(cl_data$N * mu_C1)/mean(cl_data$N * mu_C0)
    eff_ml_se_C <- sqrt(var_rr(cbind(mu_C1 - mean(mu_C1),mu_C0 - mean(mu_C0)), c(mean(mu_C1), mean(mu_C0)))/(m-5))
    eff_ml_se_I <- sqrt(var_rr(cbind(cl_data$N * (mu_C1 - mean(mu_C1))/mean(cl_data$N), cl_data$N * (mu_C0 - mean(mu_C0))/mean(cl_data$N)), 
                                c(mean(cl_data$N * mu_C1)/mean(cl_data$N), mean(cl_data$N * mu_C0)/mean(cl_data$N)))/(m-5)) 
    result_eff <- c(eff_pm_est_C, eff_pm_se_C, eff_pm_est_I, eff_pm_se_I, eff_ml_est_C, eff_ml_se_C, eff_ml_est_I, eff_ml_se_I)
    # ---------
    c(result_C, result_I, result_eff)
  }, error = function(e) { return(rep(NA, 28))})
}
tictoc::toc()

est_index_C <- c(1,3,5,7,9,21,25)
est_index_I <- c(11,13,15,17,19,23,27)
sd_index_C <- c(2,4,6,8,10,22,26)
sd_index_I <- c(12,14,16,18,20,24,28)
results[abs(results) > 1e3] <- NA
summary_results <- cbind(deltaC_bias = apply(results[est_index_C,], 1, mean, na.rm = T) - delta_C,
                         deltaC_emp_se = apply(results[est_index_C,], 1, sd, na.rm = T),
                         deltaC_avg_se = apply(results[sd_index_C,], 1, mean, na.rm = T),
                         deltaC_cp = apply(abs((results[est_index_C,]-delta_C)/results[sd_index_C,]) <= qt(0.975, m), 1, mean, na.rm = T),
                         deltaI_bias = apply(results[est_index_I,], 1, mean, na.rm = T) - delta_I,
                         deltaI_emp_se = apply(results[est_index_I,], 1, sd, na.rm = T),
                         deltaI_avg_se = apply(results[sd_index_I,], 1, mean, na.rm = T),
                         deltaI_cp = apply(abs((results[est_index_I,]-delta_I)/results[sd_index_I,]) <= qt(0.975, m), 1, mean, na.rm = T))

rownames(summary_results) <- c("unadj", "gee", "glmm", "aug-gee", "tmle", "eff-pm", "eff-ml")
round(summary_results,2)
saveRDS(results, paste0("simulation2-2-m",m,"-",simulation_scenario,".rds"))
stopCluster(cl)
