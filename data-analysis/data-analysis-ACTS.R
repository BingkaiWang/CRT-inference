rm(list=ls())
# packages
library("splines")
library(lme4)
library("boot")
library(tidyverse)
library("haven")
library(geepack) #GEE
library(lme4) # GLMM
library(CRTgeeDR) #Aug-GEE
library(SuperLearner) # TMLE
library(tmle) # TMLE
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
expit <- function(x){1/(1+exp(-x))}
se_riskdiff <- function(beta, var){
  tt <- c(exp(-beta[1]-beta[2])/(1+exp(-beta[1]-beta[2]))^2 + exp(-beta[1])/(1+exp(-beta[1]))^2, 
          exp(-beta[1]-beta[2])/(1+exp(-beta[1]-beta[2]))^2)
  sqrt(t(tt) %*% var %*% tt)
}

observed_data <- readxl::read_xlsx("data-analysis/R01_Coupon_Aim_2_Analysis_Data_MAIN_OUTCOMES_v12_20180321.xlsx") %>%
  filter(wave == 2) %>%
  dplyr::select(cluster_id = cu_code,
                A = group,
                age = patient_age_categ,
                gender = patient_gender,
                Y = malaria_test_done,
                wealth = ses_DHS_score_pooled,
                health_facility) %>%
  mutate(Y = Y ==1) %>% as.data.frame()
observed_data <- observed_data[complete.cases(observed_data),] %>% left_join(group_by(observed_data, cluster_id) %>% summarise(N=n()))
observed_data[,1] <- as.factor(observed_data[,1])
cl_data <- observed_data[,-3] %>% group_by(cluster_id) %>% summarise_all(mean)
head(observed_data)
head(cl_data)
m <- nrow(cl_data)
pi <- 0.5

# # cluster-level effects -------
# Unadjusted-Bugni
unadj_mu1 <- sum(cl_data$Y * cl_data$A * cl_data$N)/sum(cl_data$A * cl_data$N)
unadj_mu0 <- sum(cl_data$Y * (1-cl_data$A) * cl_data$N)/sum((1-cl_data$A) * cl_data$N)
unadj_est <- unadj_mu1 - unadj_mu0
unadj_se <- sd(cl_data$N/mean(cl_data$N) * (cl_data$A*(cl_data$Y-unadj_mu1)/pi - (1-cl_data$A)*(cl_data$Y - unadj_mu0)/(1-pi)))/sqrt(m-1)

# GEE
gee <- glm(Y~ . - cluster_id, weights = observed_data$N/sum(unique(observed_data$N)),
           data =observed_data, family = 'binomial')
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
gee_est <- gee_muI1_hat - gee_muI0_hat
gee_se <- sqrt(var(gee_IF[,1] - gee_IF[,2])/(m-4))
print(gee_est)
print(gee_se)

# glmm
glmm <- lmer(Y~(1|cluster_id) + . - cluster_id, data = observed_data)
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
glmm_est <- glmm_muI1_hat - glmm_muI0_hat
glmm_se <- sqrt(var(glmm_IF[,1]-glmm_IF[,2])/(m-4))


set.seed(123)
cl_tmle <- tmle(Y = cl_data$Y, A = cl_data$A, W = cl_data[,-c(1,2,4)],
                Q.SL.library = SLmethods, g1W = rep(pi, m), family = "gaussian")
cl_est <- cl_tmle$estimates$ATE$psi
cl_se <- sqrt(cl_tmle$estimates$ATE$var.psi) * sqrt(m)/sqrt(m-4)
print(cl_est)
print(cl_se)
# 
# SAAE-PM
kappa1 <- glm(A~N, family = binomial, data = cl_data) %>% predict(type = "response")
kappa0 <- 1 - kappa1
zeta1 <- predict(lm(Y~N, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
zeta0 <- predict(lm(Y~N, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
eta <- observed_data %>%
  mutate(eta1= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==1,],-c('A', 'cluster_id')), family = binomial), newdata = observed_data, type = "response")) %>%
  mutate(eta0= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==0,],-c('A', 'cluster_id')), family = binomial), newdata = observed_data, type = "response")) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_pwm_est_C <- mean(mu_C1) - mean(mu_C0)
saae_pwm_se_C <- sd(mu_C1 - mu_C0)/sqrt(m-4)
saae_pwm_est_C
saae_pwm_se_C 
# # SAAE-ML
set.seed(123)
kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("N", "health_facility")], family = "binomial",
                       SL.library = SLmethods, cvControl = list(V=5))$SL.predict
kappa0 <- 1 - kappa1
zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N", "health_facility")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N", "health_facility")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N", "health_facility")])$pred
zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N", "health_facility")])$pred
eta_fit1 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==1]), X = dplyr::select(observed_data[observed_data$A==1,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==1], SL.library = SLmethods, cvControl = list(V=2))
eta_fit0 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==0]), X = dplyr::select(observed_data[observed_data$A==0,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==0], SL.library = SLmethods, cvControl = list(V=2))
eta <- observed_data %>%
  mutate(eta1= predict(eta_fit1, newdata = dplyr::select(observed_data, -c("Y", "cluster_id","A")))$pred) %>%
  mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(observed_data, -c("Y", "cluster_id",'A')))$pred) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_ml_est_C <- mean(mu_C1) - mean(mu_C0)
saae_ml_se_C <- sd(mu_C1 - mu_C0)/sqrt(m-4)

summary_results <- matrix(data = c(unadj_est, gee_est, glmm_est, NA, cl_est, saae_pwm_est_C, saae_ml_est_C,
                                   unadj_se, gee_se, glmm_se, NA, cl_se, saae_pwm_se_C, saae_ml_se_C), nrow = 7, ncol = 2,
                          dimnames = list(c("unadj", "GEE", "LMM", "Aug-GEE", "TMLE", "PM","ML"),
                                          c("est", "sd"))) %>% as.data.frame %>%
  mutate(ci.lower = est + sd * qt(0.025, m), ci.upper = est + sd * qt(0.975, m), rov = sd^2/sd[1]^2)
round(summary_results,3)
xtable::xtable(summary_results[,-2])

# individual-level effect ----
# unadj
unadj_mu1 <- sum(cl_data$Y * cl_data$A * cl_data$N)/sum(cl_data$A * cl_data$N)
unadj_mu0 <- sum(cl_data$Y * (1-cl_data$A) * cl_data$N)/sum((1-cl_data$A) * cl_data$N)
unadj_est <- unadj_mu1 - unadj_mu0
unadj_se <- sd(cl_data$N/mean(cl_data$N) * (cl_data$A*(cl_data$Y-unadj_mu1)/pi - (1-cl_data$A)*(cl_data$Y - unadj_mu0)/(1-pi)))/sqrt(m-1)
unadj_est
unadj_se

# gee
tictoc::tic()
gee <- glm(Y~ . - cluster_id, weights = observed_data$N/sum(unique(observed_data$N)),
           data =observed_data, family = 'binomial')
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
gee_est <- gee_muI1_hat - gee_muI0_hat
gee_se <- sqrt(var(gee_IF[,1] - gee_IF[,2])/(m-4))
print(gee_est)
print(gee_se)
tictoc::toc()

# glmm
glmm <- lmer(Y~(1|cluster_id) + . - cluster_id, data = observed_data)
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
glmm_est <- glmm_muI1_hat - glmm_muI0_hat
glmm_se <- sqrt(var(glmm_IF[,1]-glmm_IF[,2])/(m-4))


#TMLE
set.seed(123)
eta_fit1 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==1]), X = dplyr::select(observed_data[observed_data$A==1,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==1], SL.library = SLmethods, cvControl = list(V=5))
eta_fit0 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==0]), X = dplyr::select(observed_data[observed_data$A==0,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==0], SL.library = SLmethods, cvControl = list(V=5))
eta <- observed_data %>%
  mutate(eta1= predict(eta_fit1, newdata = dplyr::select(observed_data, -c("Y", "cluster_id","A")))$pred) %>%
  mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(observed_data, -c("Y", "cluster_id",'A')))$pred) 
tmle_est <- mean(eta$eta1 -eta$eta0)
tmle_se <- eta %>% mutate(IC = A/pi*(Y-eta1) + eta1 - (1-A)/(1-pi)*(Y-eta0) - eta0 - tmle_est) %>%
  group_by(cluster_id) %>% summarise(IC = sum(IC))  %>% dplyr::select(IC) %>% var %>% sqrt(.)/ nrow(observed_data) * sqrt(m-4)
tmle_est
tmle_se

# SAAE-PM
kappa1 <- glm(A~N+health_facility, family = binomial, data = cl_data) %>% predict(type = "response")
kappa0 <- 1 - kappa1
zeta1 <- predict(lm(Y~N+health_facility, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
zeta0 <- predict(lm(Y~N+health_facility, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
eta <- observed_data %>%
  mutate(eta1= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==1,],-c('A', 'cluster_id')), family = binomial), newdata = observed_data, type = "response")) %>%
  mutate(eta0= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==0,],-c('A', 'cluster_id')), family = binomial), newdata = observed_data, type = "response")) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_pwm_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
saae_pwm_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - saae_pwm_est_I)/mean(cl_data$N))/(m-4))
print(saae_pwm_est_I)
print(saae_pwm_se_I)
# # SAAE-ML
set.seed(123)
kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("N", "health_facility")], family = "binomial",
                       SL.library = SLmethods, cvControl = list(V=5))$SL.predict
kappa0 <- 1 - kappa1
zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N", "health_facility")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=3))
zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N", "health_facility")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=3))
zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N", "health_facility")])$pred
zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N", "health_facility")])$pred
eta_fit1 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==1]), X = dplyr::select(observed_data[observed_data$A==1,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==1], SL.library = SLmethods, cvControl = list(V=5))
eta_fit0 <- SuperLearner(as.numeric(observed_data$Y[observed_data$A==0]), X = dplyr::select(observed_data[observed_data$A==0,], -c("Y", "cluster_id","A")), family = "binomial",
                         id = observed_data$cluster_id[observed_data$A==0], SL.library = SLmethods, cvControl = list(V=5))
eta <- observed_data %>%
  mutate(eta1= predict(eta_fit1, newdata = dplyr::select(observed_data, -c("Y", "cluster_id","A")))$pred) %>%
  mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(observed_data, -c("Y", "cluster_id",'A')))$pred) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_ml_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
saae_ml_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - saae_ml_est_I)/mean(cl_data$N))/(m-4))

summary_results <- matrix(data = c(unadj_est, gee_est, glmm_est, NA, tmle_est, saae_pwm_est_I, saae_ml_est_I,
                                   unadj_se, gee_se, glmm_se, NA, tmle_se, saae_pwm_se_I, saae_ml_se_I), nrow = 7, ncol = 2,
                          dimnames = list(c("unadj", "GEE", "LMM", "Aug-GEE", "TMLE", "PM","ML"),
                                          c("est", "sd"))) %>% as.data.frame %>%
  mutate(ci.lower = est + sd * qt(0.025, m), ci.upper = est + sd * qt(0.975, m), rov = sd^2/sd[1]^2)
round(summary_results,2)
xtable::xtable(summary_results[,-2])
